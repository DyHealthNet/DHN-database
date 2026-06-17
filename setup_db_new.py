"""
Populate the new database structure:
  - nodes               : all nodes regardless of biological type, read from a single file
                          (includes data_type and an optional grouping column)
  - edges_parametric    : association scores from a parametric test run
  - edges_nonparametric : association scores from a nonparametric test run

At least one of PARAMETRIC_EDGES_PATH / NONPARAMETRIC_EDGES_PATH must be set; if both
are set, both tables are populated independently. Queries always use one table or the
other — they are never mixed.

Nodes file columns are configured via NODES_LABEL_COLUMN / NODES_TYPE_COLUMN
(required) and NODES_DP_NAME_COLUMN / NODES_DESCRIPTION_COLUMN / NODES_XREF_COLUMN /
NODES_GROUP_COLUMN (optional) in the environment, see utils.settings.NODES_COLUMNS.

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
                            PARAMETRIC_EDGES_PATH, NONPARAMETRIC_EDGES_PATH, NODES_META_PATH,
                            NODES_COLUMNS, CHUNK_SIZE, get_logger)

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


def create_tables(engine) -> None:
    """Create nodes, edges_parametric, and edges_nonparametric if they do not exist."""
    Base.metadata.create_all(
        engine,
        tables=[Node.__table__, EdgeParametric.__table__, EdgeNonparametric.__table__],
        checkfirst=True,
    )
    logger.info("Tables created (or already existed)")


def populate_nodes(session: Session, nodes_path: str) -> None:
    """
    Read the combined nodes file and insert every node into nodes.
    Nodes already present are skipped.

    Required columns (NODES_LABEL_COLUMN, NODES_TYPE_COLUMN): unique node id, data type.
    Optional columns (NODES_DP_NAME_COLUMN, NODES_DESCRIPTION_COLUMN, NODES_XREF_COLUMN,
    NODES_GROUP_COLUMN): display name, description, xrefs, group. Missing optional
    columns are stored as NULL, and display name falls back to the node id.
    """
    start = time.perf_counter()
    logger.info(f"Loading nodes from {nodes_path}")
    df = read_data_cached(nodes_path)

    id_col = NODES_COLUMNS['unique_id']
    type_col = NODES_COLUMNS['data_type']
    name_col = NODES_COLUMNS['display_name']
    desc_col = NODES_COLUMNS['description']
    xref_col = NODES_COLUMNS['xref']
    group_col = NODES_COLUMNS['group']

    for col, env_var in [(id_col, 'NODES_LABEL_COLUMN'), (type_col, 'NODES_TYPE_COLUMN')]:
        if not col:
            raise ValueError(f"{env_var} must be set in the environment")
        if col not in df.columns:
            raise ValueError(f"Column '{col}' (from {env_var}) not found in {nodes_path}")

    df = df.dropna(subset=[id_col])

    nodes = pd.DataFrame({
        'node_id': df[id_col],
        'display_name': df[name_col] if name_col and name_col in df.columns else df[id_col],
        'data_type': df[type_col],
        'node_group': df[group_col] if group_col and group_col in df.columns else None,
        'description': df[desc_col] if desc_col and desc_col in df.columns else None,
        'xrefs': df[xref_col] if xref_col and xref_col in df.columns else None,
    })

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
    populate_nodes(db_session, NODES_META_PATH)

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
