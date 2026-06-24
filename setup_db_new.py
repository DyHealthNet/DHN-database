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
display_name/description/xref/group per node where available. These are optional,
per-source columns: DATA_DP_NAME_COLUMNS / DATA_DESCRIPTION_COLUMNS /
DATA_XREF_COLUMNS / DATA_GROUP_COLUMNS. If set, each must have the same number of
comma-separated entries as DATA_META_PATHS (use an empty entry to skip a source that
doesn't have that column). Left unset entirely, nodes just get no enrichment
(display_name falls back to the node id, the rest stay NULL) - the program never fails
because these are missing.

Expected scores CSV format (produced by the association-score package):
    index, label1, label2, raw-P, raw-E, test_type

CSV/TSV inputs are cached as parquet files alongside the source file for faster
repeated reads (mirrors DHN-backend's startup_utils.check_files_and_return). The
scores file is streamed from the parquet cache in chunks of CHUNK_SIZE rows, so
arbitrarily large files can be inserted without loading them fully into memory.

Run directly:
    python setup_db_new.py
"""
import os
import time
import pandas as pd
import pyarrow.parquet as pq
from io import StringIO

from sqlalchemy import URL, create_engine, text
from sqlalchemy.orm import sessionmaker, Session

from utils.models import Base, Node, EdgeParametric, EdgeNonparametric
from utils.settings import (DB_USER, DB_PASSWORD, DB_HOST, DB_PORT, DB_NAME,
                            PARAMETRIC_EDGES_PATH, NONPARAMETRIC_EDGES_PATH, CHUNK_SIZE, get_logger,
                            DATA_ROOT, DATA_META_PATHS, DATA_LABEL_COLUMNS, DATA_TYPE_COLUMNS,
                            DATA_DP_NAME_COLUMNS, DATA_DESCRIPTION_COLUMNS,
                            DATA_XREF_COLUMNS, DATA_GROUP_COLUMNS)

logger = get_logger(__name__)

SCORES_COLUMNS = ['label1', 'label2', 'raw-P', 'raw-E', 'test_type']

_VALID_EDGE_TABLES = {'edges_parametric', 'edges_nonparametric'}


def _ensure_parquet_cache(path: str) -> str:
    """
    Return the path to a parquet version of the given CSV/TSV/Parquet file,
    creating the cache (alongside the source file) if it doesn't exist or is stale.
    """
    file_name, ending = os.path.splitext(path)
    ending = ending.lower()

    if ending == '.parquet':
        return path
    if ending not in ['.csv', '.tsv']:
        raise ValueError(f"Unsupported file format: {ending}. Only CSV, TSV and Parquet files are supported.")

    parquet_path = f"{file_name}.parquet"
    if not os.path.exists(parquet_path) or os.path.getmtime(path) > os.path.getmtime(parquet_path):
        logger.debug(f"Creating parquet cache for {path}")
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
    """Bulk-insert a DataFrame into a table via PostgreSQL COPY."""
    buffer = StringIO()
    df[columns].to_csv(buffer, index=False, header=False, na_rep='')
    buffer.seek(0)

    raw_conn = session.connection().connection
    with raw_conn.cursor() as cursor:
        col_list = ', '.join(columns)
        copy_sql = f"COPY {table} ({col_list}) FROM STDIN WITH (FORMAT CSV, NULL '')"
        cursor.copy_expert(copy_sql, buffer)
    raw_conn.commit()


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


def _load_node_source(meta_path: str, label_column: str, type_column: str,
                      name_column: str, desc_column: str, xref_column: str, group_column: str) -> pd.DataFrame:
    """
    Read one meta file into a [node_id, display_name, data_type, node_group,
    description, xrefs] frame. name/desc/xref/group_column are optional (pass "" to
    skip); display_name falls back to node_id and the rest stay NULL when skipped or
    not present in this file.
    """
    df = read_data_cached(meta_path)
    if label_column not in df.columns:
        raise ValueError(f"Column '{label_column}' not found in {meta_path}")

    df = df.rename(columns={label_column: 'node_id'})
    df = _apply_type_column(df, type_column)
    df = df.dropna(subset=['node_id'])

    return pd.DataFrame({
        'node_id': df['node_id'],
        'display_name': df[name_column] if name_column and name_column in df.columns else df['node_id'],
        'data_type': df['data_type'],
        'node_group': df[group_column] if group_column and group_column in df.columns else None,
        'description': df[desc_column] if desc_column and desc_column in df.columns else None,
        'xrefs': df[xref_column] if xref_column and xref_column in df.columns else None,
    })


def build_combined_nodes() -> pd.DataFrame:
    """
    Build the combined nodes table from DATA_META_PATHS, mirroring DHN-backend's
    network.utils.data_manager.combine_data(): one entry per data source, each
    contributing label/type (required). DATA_DP_NAME_COLUMNS/DATA_DESCRIPTION_COLUMNS/
    DATA_XREF_COLUMNS/DATA_GROUP_COLUMNS optionally add display_name/description/xref/
    group per source - the backend doesn't need these for scoring, but the DB wants
    them where available. Missing optional columns are stored as NULL; display name
    falls back to the node id.
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

    frames = []
    for meta_path, label_column, type_column, name_column, desc_column, xref_column, group_column in zip(
        meta_paths, label_columns, type_columns, name_columns, desc_columns, xref_columns, group_columns
    ):
        meta_path = _resolve_path(meta_path, DATA_ROOT)
        frames.append(_load_node_source(
            meta_path, label_column, type_column, name_column, desc_column, xref_column, group_column
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
                    ['node_id', 'display_name', 'data_type', 'node_group', 'description', 'xrefs'])
    logger.info(f"Inserted {len(nodes)} nodes in {time.perf_counter() - start:.2f}s")


def insert_scores(session: Session, edges_path: str, table: str,
                  truncate: bool = False, chunk_size: int = CHUNK_SIZE) -> None:
    """
    Read a scores file and insert into the given edge table (edges_parametric or
    edges_nonparametric), processing it in chunks so arbitrarily large files don't
    need to fit in memory.

    Expected columns: label1, label2, raw-P, raw-E, test_type.
    Rows missing raw-P or whose node IDs are not present in nodes are skipped.

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

    total_inserted = 0
    total_skipped = 0
    for df in iter_data_chunks(edges_path, chunk_size):
        df = df.loc[:, [c for c in df.columns if not c.startswith('Unnamed')]]

        missing_cols = [c for c in SCORES_COLUMNS if c not in df.columns]
        if missing_cols:
            raise ValueError(f"Scores file is missing expected columns: {missing_cols}")

        df = df[SCORES_COLUMNS].dropna(subset=['label1', 'label2', 'raw-P'])

        fk_mask = df['label1'].isin(known_nodes) & df['label2'].isin(known_nodes)
        total_skipped += int((~fk_mask).sum())
        df = df[fk_mask]

        if df.empty:
            continue

        df = df.rename(columns={'label1': 'node_id_1', 'label2': 'node_id_2', 'raw-P': 'p_value', 'raw-E': 'effect_size'})

        _copy_dataframe(session, df, table, ['node_id_1', 'node_id_2', 'p_value', 'effect_size', 'test_type'])
        total_inserted += len(df)

    elapsed = time.perf_counter() - start
    if total_skipped:
        logger.warning(f"Skipped {total_skipped} edges with node IDs not in nodes")
    if total_inserted:
        logger.info(f"Inserted {total_inserted} rows into {table} in {elapsed:.2f}s")
    else:
        logger.warning(f"No valid edges to insert into {table} ({elapsed:.2f}s)")


if __name__ == '__main__':
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

    total_start = time.perf_counter()
    create_tables(engine)
    populate_nodes(db_session)

    if PARAMETRIC_EDGES_PATH:
        insert_scores(db_session, PARAMETRIC_EDGES_PATH, table='edges_parametric')
    else:
        logger.info("PARAMETRIC_EDGES_PATH not set, skipping edges_parametric")

    if NONPARAMETRIC_EDGES_PATH:
        insert_scores(db_session, NONPARAMETRIC_EDGES_PATH, table='edges_nonparametric')
    else:
        logger.info("NONPARAMETRIC_EDGES_PATH not set, skipping edges_nonparametric")

    db_session.close()
    logger.info(f"Database population complete in {time.perf_counter() - total_start:.2f}s.")
