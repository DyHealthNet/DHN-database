"""
Populate the new database structure:
  - nodes               : all nodes regardless of biological type, combined from any
                          number of data sources (includes data_type and an optional
                          grouping column)
  - edges_parametric    : association scores from a parametric test run
  - edges_nonparametric : association scores from a nonparametric test run

At least one of PARAMETRIC_EDGES_PATH / NONPARAMETRIC_EDGES_PATH must be set; if both
are set, both tables are populated independently. Queries always use one table or the
other — they are never mixed.

Nodes are combined from any number of data sources (e.g. phenotypes, proteins,
metabolites), mirroring DHN-backend's network.utils.data_manager.combine_data():
configured via the comma-separated DATA_META_PATHS / DATA_LABEL_COLUMNS /
DATA_TYPE_COLUMNS, which must all have the same number of entries (one per source, at
least one required). DATA_ROOT, if set, is prepended to relative entries in
DATA_META_PATHS. Each source contributes label/type (required) to the nodes table,
matching what the backend uses for scoring.

The backend doesn't need them for scoring, but the DB additionally wants
display_name/description/xref/group/subgroup per node where available. These are
optional, per-source columns: DATA_DP_NAME_COLUMNS / DATA_DESCRIPTION_COLUMNS /
DATA_XREF_COLUMNS / DATA_GROUP_COLUMNS / DATA_SUBGROUP_COLUMNS. If set, each must have
the same number of comma-separated entries as DATA_META_PATHS (use an empty entry to
skip a source that doesn't have that column). Left unset entirely, nodes just get no
enrichment (display_name falls back to the node id, the rest stay NULL) - the program
never fails because these are missing.

Expected scores CSV format (produced by the association-score package):
    index, label1, label2, raw-P, raw-E, test_type

CSV/TSV inputs are cached as parquet files alongside the source file for faster
repeated reads (mirrors DHN-backend's startup_utils.check_files_and_return). The
scores file is streamed from the parquet cache in chunks of CHUNK_SIZE rows, so
arbitrarily large files can be inserted without loading them fully into memory.

Run directly:
    python setup_db_new.py
"""
import argparse
import os
import time
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pc
import pyarrow.parquet as pq
from io import BytesIO

from sqlalchemy import URL, create_engine, text
from sqlalchemy.orm import sessionmaker, Session

from utils.models import Base, Node, EdgeParametric, EdgeNonparametric
from utils.settings import (DB_USER, DB_PASSWORD, DB_HOST, DB_PORT, DB_NAME,
                            PARAMETRIC_EDGES_PATH, NONPARAMETRIC_EDGES_PATH, CHUNK_SIZE, get_logger,
                            DATA_ROOT, DATA_META_PATHS, DATA_LABEL_COLUMNS, DATA_TYPE_COLUMNS,
                            DATA_DP_NAME_COLUMNS, DATA_DESCRIPTION_COLUMNS,
                            DATA_XREF_COLUMNS, DATA_GROUP_COLUMNS, DATA_SUBGROUP_COLUMNS, DATA_DIR)
from utils.database_views import add_indexes_new

logger = get_logger(__name__)

SCORES_COLUMNS = ['label1', 'label2', 'raw-P', 'raw-E', 'test_type']

_VALID_EDGE_TABLES = {'edges_parametric', 'edges_nonparametric'}


def _ensure_parquet_cache(path: str) -> str:
    """
    Return the path to a parquet version of the given CSV/TSV/Parquet file,
    creating the cache if it doesn't exist or is stale.

    Cached under DATA_DIR (a directory this process is expected to own) rather than
    alongside the source file, since source files often live on shared, read-only
    mounts (e.g. DATA_ROOT pointing at another user's data) this process has no
    write access to. The cache filename is just the source's basename, so two
    sources with the same basename in different directories will collide.
    Falls back to caching alongside the source file if DATA_DIR isn't configured.
    """
    file_name, ending = os.path.splitext(path)
    ending = ending.lower()

    if ending == '.parquet':
        return path
    if ending not in ['.csv', '.tsv']:
        raise ValueError(f"Unsupported file format: {ending}. Only CSV, TSV and Parquet files are supported.")

    if DATA_DIR:
        os.makedirs(DATA_DIR, exist_ok=True)
        cache_name = os.path.basename(file_name) + '.parquet'
        parquet_path = os.path.join(DATA_DIR, cache_name)
    else:
        parquet_path = f"{file_name}.parquet"

    if not os.path.exists(parquet_path) or os.path.getmtime(path) > os.path.getmtime(parquet_path):
        logger.debug(f"Creating parquet cache for {path} at {parquet_path}")
        sep = ',' if ending == '.csv' else '\t'
        df = pd.read_csv(path, sep=sep, low_memory=False)
        df.to_parquet(parquet_path)

    return parquet_path


def read_data_cached(path: str) -> pd.DataFrame:
    """
    Read a CSV/TSV/Parquet file in full. CSV/TSV files are cached as a parquet
    file alongside the source file; the cache is rebuilt if the source is newer.
    """
    parquet_path = _ensure_parquet_cache(path)
    logger.debug(f"Reading {parquet_path}")
    return pd.read_parquet(parquet_path)


def iter_data_chunks(path: str, chunk_size: int = CHUNK_SIZE):
    """
    Yield a CSV/TSV/Parquet file as successive DataFrame chunks of at most
    chunk_size rows each, via the parquet cache (see read_data_cached).
    """
    parquet_path = _ensure_parquet_cache(path)
    parquet_file = pq.ParquetFile(parquet_path)
    for batch in parquet_file.iter_batches(batch_size=chunk_size):
        yield batch.to_pandas()


def _copy_dataframe(session: Session, df: pd.DataFrame, table: str, columns: list[str]) -> None:
    """
    Bulk-insert a DataFrame into a table via PostgreSQL COPY. Serializes via PyArrow
    (much faster than pandas' to_csv for large frames) and disables triggers for the
    duration of the COPY to skip per-row FK validation - callers are expected to have
    already filtered rows so their FKs are valid (see insert_scores' known_nodes check).
    """
    buffer = BytesIO()
    arrow_table = pa.Table.from_pandas(df[columns], preserve_index=False)
    pc.write_csv(arrow_table, buffer, write_options=pc.WriteOptions(include_header=False, quoting_style=None))
    buffer.seek(0)

    raw_conn = session.connection().connection
    with raw_conn.cursor() as cursor:
        col_list = ', '.join(columns)
        copy_sql = f"COPY {table} ({col_list}) FROM STDIN WITH (FORMAT CSV, NULL '')"
        cursor.execute(f"ALTER TABLE {table} DISABLE TRIGGER ALL")
        cursor.copy_expert(copy_sql, buffer)
        cursor.execute(f"ALTER TABLE {table} ENABLE TRIGGER ALL")
    raw_conn.commit()


def _copy_edges_ignoring_conflicts(raw_conn, df: pd.DataFrame, columns: list[str],
                                    table: str, staging_table: str) -> int:
    """
    Load `df` into `table` (edges_parametric/edges_nonparametric) via a
    session-scoped staging table plus a single set-based `INSERT ... ON CONFLICT
    DO NOTHING`, canonicalizing node_id_1/node_id_2 order (LEAST/GREATEST) as
    part of that same INSERT.

    Canonicalizing here in SQL -- not earlier in pandas -- matters: Postgres's
    default locale collation doesn't sort strings the same way Python's plain
    Unicode code-point comparison does (verified empirically against this
    database: 'hcu' < 'M01A' under its en_US.utf8 collation, but 'M01A' < 'hcu'
    in Python). The table's UniqueConstraint on (node_id_1, node_id_2) is itself
    collation-based, so only a same-collation LEAST/GREATEST is guaranteed to
    agree with it -- a Python-side min/max can canonicalize a pair to the
    opposite order from what's already stored, silently defeating the
    constraint for exactly the mixed-case ids where the two orderings disagree.

    Returns the number of rows actually inserted.
    """
    buffer = BytesIO()
    arrow_table = pa.Table.from_pandas(df[columns], preserve_index=False)
    pc.write_csv(arrow_table, buffer, write_options=pc.WriteOptions(include_header=False, quoting_style=None))
    buffer.seek(0)

    col_list = ', '.join(columns)
    with raw_conn.cursor() as cursor:
        cursor.execute(f"TRUNCATE {staging_table}")
        cursor.copy_expert(f"COPY {staging_table} ({col_list}) FROM STDIN WITH (FORMAT CSV, NULL '')", buffer)
        # Unique-constraint enforcement (what ON CONFLICT below relies on) is
        # index-based, not trigger-based, so disabling triggers here still only
        # skips the (already pre-validated via known_nodes) FK checks, same as
        # the plain _copy_dataframe path.
        cursor.execute(f"ALTER TABLE {table} DISABLE TRIGGER ALL")
        cursor.execute(f"""
            INSERT INTO {table} (node_id_1, node_id_2, p_value, effect_size, test_type)
            SELECT LEAST(node_id_1, node_id_2), GREATEST(node_id_1, node_id_2),
                   p_value, effect_size, test_type
            FROM {staging_table}
            ON CONFLICT (node_id_1, node_id_2) DO NOTHING
        """)
        inserted = cursor.rowcount
        cursor.execute(f"ALTER TABLE {table} ENABLE TRIGGER ALL")
    raw_conn.commit()
    return inserted


def _parse_list_env(value: str | None) -> list[str]:
    """Parse a comma-separated env var into a list of stripped, non-empty strings."""
    return [item.strip() for item in (value or "").split(",") if item.strip()]


def _parse_aligned_list_env(value: str | None, expected_length: int, env_var: str) -> list[str]:
    """
    Parse an optional comma-separated env var into a list aligned 1:1 with the data
    sources. Unlike _parse_list_env, empty entries are kept (not dropped) so position
    still lines up with DATA_META_PATHS - use an empty entry to skip a source that
    doesn't have that column. If unset entirely, returns `expected_length` empty
    strings (no enrichment for any source).
    """
    if not value:
        return [""] * expected_length

    parts = [item.strip() for item in value.split(",")]
    if len(parts) != expected_length:
        raise ValueError(
            f"{env_var} has {len(parts)} comma-separated entries, expected {expected_length} "
            f"(one per DATA_META_PATHS entry; use an empty entry to skip a source)."
        )
    return parts


def _resolve_path(path: str, root: str | None) -> str:
    """Join a relative path with `root` (if set); absolute paths are returned unchanged."""
    if root and not os.path.isabs(path):
        return os.path.join(root, path)
    return path


def _apply_type_column(meta: pd.DataFrame, type_column: str) -> pd.DataFrame:
    """
    Rename `type_column` to 'data_type' if it is an existing column in `meta`.
    Otherwise treat `type_column` as a literal type value applied to every row.
    Mirrors DHN-backend's network.utils.data_manager._apply_type_column.
    """
    if type_column in meta.columns:
        return meta.rename(columns={type_column: 'data_type'})
    meta = meta.copy()
    meta['data_type'] = type_column
    return meta


def _apply_literal_or_column(meta: pd.DataFrame, column_or_literal: str) -> pd.Series | None:
    """
    Return `column_or_literal`'s values if it names an existing column in `meta`;
    otherwise treat it as a literal value applied to every row (e.g. a fixed group
    name for a source with no per-row grouping of its own). Returns None if
    `column_or_literal` is falsy (the attribute isn't configured for this source).
    Mirrors DHN-backend's network.utils.data_manager._apply_literal_or_column.
    """
    if not column_or_literal:
        return None
    if column_or_literal in meta.columns:
        return meta[column_or_literal]
    return pd.Series(column_or_literal, index=meta.index)


def _load_node_source(meta_path: str, label_column: str, type_column: str,
                      name_column: str, desc_column: str, xref_column: str, group_column: str,
                      subgroup_column: str) -> pd.DataFrame:
    """
    Read one meta file into a [node_id, display_name, data_type, node_group,
    node_subgroup, description, xrefs] frame. name/desc/xref/group/subgroup_column are
    optional (pass "" to skip) and, like type_column, may each be either a column name
    in this meta file or a literal value applied to every row. display_name falls back
    to node_id and the rest stay NULL when skipped.
    """
    df = read_data_cached(meta_path)
    if label_column not in df.columns:
        raise ValueError(f"Column '{label_column}' not found in {meta_path}")

    df = df.rename(columns={label_column: 'node_id'})
    df = _apply_type_column(df, type_column)
    df = df.dropna(subset=['node_id'])

    display_name = _apply_literal_or_column(df, name_column)

    return pd.DataFrame({
        'node_id': df['node_id'],
        'display_name': display_name if display_name is not None else df['node_id'],
        'data_type': df['data_type'],
        'node_group': _apply_literal_or_column(df, group_column),
        'node_subgroup': _apply_literal_or_column(df, subgroup_column),
        'description': _apply_literal_or_column(df, desc_column),
        'xrefs': _apply_literal_or_column(df, xref_column),
    })


def build_combined_nodes() -> pd.DataFrame:
    """
    Build the combined nodes table from DATA_META_PATHS, mirroring DHN-backend's
    network.utils.data_manager.combine_data(): one entry per data source, each
    contributing label/type (required). DATA_DP_NAME_COLUMNS/DATA_DESCRIPTION_COLUMNS/
    DATA_XREF_COLUMNS/DATA_GROUP_COLUMNS/DATA_SUBGROUP_COLUMNS optionally add
    display_name/description/xref/group/subgroup per source - the backend doesn't need
    these for scoring, but the DB wants them where available. Missing optional columns
    are stored as NULL; display name falls back to the node id.
    """
    meta_paths = _parse_list_env(DATA_META_PATHS)
    label_columns = _parse_list_env(DATA_LABEL_COLUMNS)
    type_columns = _parse_list_env(DATA_TYPE_COLUMNS)

    if not (len(meta_paths) == len(label_columns) == len(type_columns)):
        raise ValueError(
            "DATA_META_PATHS, DATA_LABEL_COLUMNS and DATA_TYPE_COLUMNS must all have "
            "the same number of comma-separated entries."
        )
    if not meta_paths:
        raise ValueError(
            "No node source configured. Set DATA_META_PATHS, DATA_LABEL_COLUMNS and DATA_TYPE_COLUMNS."
        )

    name_columns = _parse_aligned_list_env(DATA_DP_NAME_COLUMNS, len(meta_paths), 'DATA_DP_NAME_COLUMNS')
    desc_columns = _parse_aligned_list_env(DATA_DESCRIPTION_COLUMNS, len(meta_paths), 'DATA_DESCRIPTION_COLUMNS')
    xref_columns = _parse_aligned_list_env(DATA_XREF_COLUMNS, len(meta_paths), 'DATA_XREF_COLUMNS')
    group_columns = _parse_aligned_list_env(DATA_GROUP_COLUMNS, len(meta_paths), 'DATA_GROUP_COLUMNS')
    subgroup_columns = _parse_aligned_list_env(DATA_SUBGROUP_COLUMNS, len(meta_paths), 'DATA_SUBGROUP_COLUMNS')

    frames = []
    for meta_path, label_column, type_column, name_column, desc_column, xref_column, group_column, subgroup_column in zip(
        meta_paths, label_columns, type_columns, name_columns, desc_columns, xref_columns, group_columns, subgroup_columns
    ):
        meta_path = _resolve_path(meta_path, DATA_ROOT)
        frames.append(_load_node_source(
            meta_path, label_column, type_column, name_column, desc_column, xref_column, group_column, subgroup_column
        ))

    combined = pd.concat(frames, ignore_index=True)
    duplicates = combined['node_id'].duplicated()
    if duplicates.any():
        logger.warning(f"Dropping {duplicates.sum()} duplicate node id(s) found across node sources.")
        combined = combined[~duplicates]

    return combined


def create_tables(engine) -> None:
    """Create nodes, edges_parametric, and edges_nonparametric if they do not exist."""
    Base.metadata.create_all(
        engine,
        tables=[Node.__table__, EdgeParametric.__table__, EdgeNonparametric.__table__],
        checkfirst=True,
    )
    logger.info("Tables created (or already existed)")


def drop_tables(engine) -> None:
    """
    Drop nodes, edges_parametric, and edges_nonparametric (CASCADE, so anything
    referencing them, e.g. their indexes/FKs, goes with them).

    Also drops every per-context table (edges_parametric_<context_id> /
    edges_nonparametric_<context_id>, created on demand by DHN-backend's
    network/contexts/contexts.py when a user builds a context) since those were
    computed against the data being replaced and would otherwise silently go stale.
    context and user_context (DHN-backend's network.models.Context /
    UserContextLink tables) are cleared -- not dropped, the Django app still owns
    their schema -- since every context they reference no longer exists.

    Only ever called when --overwrite is passed; create_tables() recreates
    nodes/edges_parametric/edges_nonparametric fresh right after, and
    add_indexes_new() rebuilds their indexes at the end of the run.
    """
    logger.warning(
        "--overwrite passed: dropping edges_parametric, edges_nonparametric, nodes, and all "
        "per-context edges_*_<context_id> tables; clearing context and user_context."
    )
    with engine.connect() as conn:
        context_tables = conn.execute(text("""
            SELECT table_name FROM information_schema.tables
            WHERE table_name ~ '^edges_(parametric|nonparametric)_[0-9]+$'
        """)).scalars().all()
        for table in context_tables:
            conn.execute(text(f"DROP TABLE IF EXISTS {table} CASCADE"))
        if context_tables:
            logger.info(f"Dropped {len(context_tables)} per-context edges table(s): {context_tables}")

        existing_link_tables = set(conn.execute(text("""
            SELECT table_name FROM information_schema.tables
            WHERE table_name IN ('user_context', 'context')
        """)).scalars().all())
        for table in ('user_context', 'context'):
            if table in existing_link_tables:
                conn.execute(text(f"DELETE FROM {table}"))
            else:
                logger.debug(f"{table} table not present, skipping clear")

        for table in ('edges_parametric', 'edges_nonparametric', 'nodes'):
            conn.execute(text(f"DROP TABLE IF EXISTS {table} CASCADE"))
        conn.commit()


def populate_nodes(session: Session) -> None:
    """
    Build the combined nodes table (see build_combined_nodes) and insert every node into
    nodes. Nodes already present are skipped.
    """
    start = time.perf_counter()
    nodes = build_combined_nodes()

    existing_ids = {row[0] for row in session.query(Node.node_id).all()}
    nodes = nodes[~nodes['node_id'].isin(existing_ids)]

    if nodes.empty:
        logger.info(f"No new nodes to insert ({time.perf_counter() - start:.2f}s)")
        return

    _copy_dataframe(session, nodes, 'nodes',
                    ['node_id', 'display_name', 'data_type', 'node_group', 'node_subgroup', 'description', 'xrefs'])
    logger.info(f"Inserted {len(nodes)} nodes in {time.perf_counter() - start:.2f}s")


def insert_scores(session: Session, edges_path: str, table: str,
                  truncate: bool = False, chunk_size: int = CHUNK_SIZE) -> None:
    """
    Read a scores file and insert into the given edge table (edges_parametric or
    edges_nonparametric), processing it in chunks so arbitrarily large files don't
    need to fit in memory.

    Expected columns: label1, label2, raw-P, raw-E, test_type.
    Rows missing raw-P or whose node IDs are not present in nodes are skipped.

    node_id_1/node_id_2 are canonicalized (smaller node_id first) so a pair already
    loaded in either order is recognized as the same undirected edge: rows that
    collide with an already-present pair (via the table's UniqueConstraint, see
    edges_parametric_pair_unique/edges_nonparametric_pair_unique in utils/models.py)
    are skipped via `ON CONFLICT DO NOTHING` rather than inserted as a duplicate.
    This makes re-running against an already-populated table a cheap no-op instead
    of appending a second full copy of every edge -- the dedup check itself is a
    Postgres index lookup done set-at-a-time per chunk, not a Python-side membership
    test against the whole existing table, so it stays cheap as the table grows.

    :param session:     SQLAlchemy session connected to the target database.
    :param edges_path:  Path to the scores CSV/TSV/Parquet file.
    :param table:       Target table name ('edges_parametric' or 'edges_nonparametric').
    :param truncate:    If True, clear the target table before inserting.
    :param chunk_size:  Number of rows to read and insert per batch.
    """
    if table not in _VALID_EDGE_TABLES:
        raise ValueError(f"Invalid edge table '{table}'. Must be one of: {_VALID_EDGE_TABLES}")

    start = time.perf_counter()
    logger.info(f"Loading scores from {edges_path} into {table} in chunks of {chunk_size:,} rows")

    if truncate:
        session.execute(text(f"TRUNCATE TABLE {table}"))
        session.commit()
        logger.info(f"Truncated {table}")

    # FK safety: only insert edges whose both node IDs exist in nodes
    known_nodes = {row[0] for row in session.query(Node.node_id).all()}

    edge_columns = ['node_id_1', 'node_id_2', 'p_value', 'effect_size', 'test_type']
    staging_table = f"_{table}_staging"
    raw_conn = session.connection().connection
    with raw_conn.cursor() as cursor:
        # No id/constraints here on purpose -- this only ever holds one chunk at a
        # time before being merged into `table` via ON CONFLICT, see
        # _copy_dataframe_ignoring_conflicts.
        cursor.execute(f"""
            CREATE TEMP TABLE IF NOT EXISTS {staging_table} (
                node_id_1 varchar, node_id_2 varchar,
                p_value double precision, effect_size double precision, test_type varchar
            ) ON COMMIT PRESERVE ROWS
        """)
    raw_conn.commit()

    total_inserted = 0
    total_duplicate = 0
    total_skipped = 0
    total_missing = 0
    for df in iter_data_chunks(edges_path, chunk_size):
        df = df.loc[:, [c for c in df.columns if not c.startswith('Unnamed')]]

        missing_cols = [c for c in SCORES_COLUMNS if c not in df.columns]
        if missing_cols:
            raise ValueError(f"Scores file is missing expected columns: {missing_cols}")

        df = df[SCORES_COLUMNS]
        na_mask = df[['label1', 'label2', 'raw-P']].isna().any(axis=1)
        total_missing += int(na_mask.sum())
        df = df[~na_mask]

        fk_mask = df['label1'].isin(known_nodes) & df['label2'].isin(known_nodes)
        total_skipped += int((~fk_mask).sum())
        df = df[fk_mask]

        if df.empty:
            continue

        df = df.rename(columns={'label1': 'node_id_1', 'label2': 'node_id_2', 'raw-P': 'p_value', 'raw-E': 'effect_size'})

        # Pair order (smaller node_id first) is canonicalized in SQL inside
        # _copy_edges_ignoring_conflicts, not here -- see its docstring for why.
        chunk_len = len(df)
        inserted = _copy_edges_ignoring_conflicts(raw_conn, df, edge_columns, table, staging_table)
        total_inserted += inserted
        total_duplicate += chunk_len - inserted

    with raw_conn.cursor() as cursor:
        cursor.execute(f"DROP TABLE IF EXISTS {staging_table}")
    raw_conn.commit()

    elapsed = time.perf_counter() - start
    if total_missing:
        logger.warning(f"Skipped {total_missing} edges with missing label1/label2/raw-P")
    if total_skipped:
        logger.warning(f"Skipped {total_skipped} edges with node IDs not in nodes")
    if total_duplicate:
        logger.info(f"Skipped {total_duplicate} edges already present in {table}")
    if total_inserted:
        logger.info(f"Inserted {total_inserted} rows into {table} in {elapsed:.2f}s")
    else:
        logger.warning(f"No new edges to insert into {table} ({elapsed:.2f}s)")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '--overwrite', action='store_true',
        help="Drop nodes, edges_parametric, edges_nonparametric, and all per-context edges_*_<id> "
             "tables before repopulating, and clear the context/user_context tables (since every "
             "context they reference becomes invalid), instead of the default behavior (existing "
             "nodes/edges/contexts are left as-is, only new rows are added)."
    )
    args = parser.parse_args()

    if not PARAMETRIC_EDGES_PATH and not NONPARAMETRIC_EDGES_PATH:
        raise EnvironmentError(
            "At least one of PARAMETRIC_EDGES_PATH or NONPARAMETRIC_EDGES_PATH must be set. "
            "The platform cannot run without any edges."
        )

    url_obj = URL.create(
        "postgresql",
        username=DB_USER,
        password=DB_PASSWORD,
        host=DB_HOST,
        port=DB_PORT,
        database=DB_NAME,
    )
    engine = create_engine(url_obj)
    SessionFactory = sessionmaker(bind=engine)
    db_session = SessionFactory()

    if args.overwrite:
        drop_tables(engine)

    total_start = time.perf_counter()
    create_tables(engine)
    populate_nodes(db_session)

    if PARAMETRIC_EDGES_PATH:
        insert_scores(db_session, PARAMETRIC_EDGES_PATH, table='edges_parametric', truncate=args.overwrite)
    else:
        logger.info("PARAMETRIC_EDGES_PATH not set, skipping edges_parametric")

    if NONPARAMETRIC_EDGES_PATH:
        insert_scores(db_session, NONPARAMETRIC_EDGES_PATH, table='edges_nonparametric', truncate=args.overwrite)
    else:
        logger.info("NONPARAMETRIC_EDGES_PATH not set, skipping edges_nonparametric")

    logger.info("Adding indexes...")
    add_indexes_new(db_session, engine)

    db_session.close()
    logger.info(f"Database population complete in {time.perf_counter() - total_start:.2f}s.")
