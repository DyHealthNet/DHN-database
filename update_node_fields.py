"""
Update display_name/xrefs for existing node rows in place, without touching node_id,
edges, or indexes -- e.g. after correcting a source's UniProt mapping post-import.
The complement of setup_db_new.py's populate_nodes (insert-only): this is
update-only, and only ever touches nodes whose node_id already exists.

Usage:
    # 1. Export current rows so you can merge fixes into xrefs without losing
    #    whatever else is already stored there.
    python update_node_fields.py --export current_proteins.csv --data-type protein

    # 2. Edit current_proteins.csv (fix display_name/xrefs), then apply it.
    python update_node_fields.py --apply fixed_proteins.csv
"""
import argparse
import csv
from io import BytesIO

import pandas as pd
import pyarrow as pa
import pyarrow.csv as pc
from sqlalchemy import URL, create_engine
from sqlalchemy.orm import sessionmaker, Session

from utils.models import Node
from utils.settings import DB_USER, DB_PASSWORD, DB_HOST, DB_PORT, DB_NAME, get_logger

logger = get_logger(__name__)


def _get_engine():
    url_obj = URL.create(
        "postgresql",
        username=DB_USER,
        password=DB_PASSWORD,
        host=DB_HOST,
        port=DB_PORT,
        database=DB_NAME,
    )
    return create_engine(url_obj)


def export_nodes(session: Session, output_path: str, data_type: str | None) -> None:
    query = session.query(Node.node_id, Node.display_name, Node.xrefs)
    if data_type:
        query = query.filter(Node.data_type == data_type)
    rows = query.all()

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['node_id', 'display_name', 'xrefs'])
        writer.writerows(rows)
    logger.info(f"Exported {len(rows)} node(s) to {output_path}")


def apply_fixes(session: Session, input_path: str) -> None:
    fixes = pd.read_csv(input_path, dtype=str, keep_default_na=False)
    missing_cols = [c for c in ('node_id', 'display_name', 'xrefs') if c not in fixes.columns]
    if missing_cols:
        raise ValueError(f"{input_path} is missing expected column(s): {missing_cols}")
    fixes = fixes[['node_id', 'display_name', 'xrefs']].replace({'': None})

    existing_ids = {row[0] for row in session.query(Node.node_id).all()}
    unknown = set(fixes['node_id']) - existing_ids
    if unknown:
        preview = sorted(unknown)[:10]
        logger.warning(f"Skipping {len(unknown)} node_id(s) not present in nodes: {preview}{'...' if len(unknown) > 10 else ''}")
        fixes = fixes[fixes['node_id'].isin(existing_ids)]

    if fixes.empty:
        logger.info("Nothing to update.")
        return

    raw_conn = session.connection().connection
    staging_table = '_node_fixes_staging'
    buffer = BytesIO()
    arrow_table = pa.Table.from_pandas(fixes, preserve_index=False)
    pc.write_csv(arrow_table, buffer, write_options=pc.WriteOptions(include_header=False, quoting_style=None))
    buffer.seek(0)

    with raw_conn.cursor() as cursor:
        cursor.execute(f"""
            CREATE TEMP TABLE IF NOT EXISTS {staging_table} (
                node_id varchar, display_name varchar, xrefs varchar
            ) ON COMMIT PRESERVE ROWS
        """)
        cursor.copy_expert(
            f"COPY {staging_table} (node_id, display_name, xrefs) FROM STDIN WITH (FORMAT CSV, NULL '')", buffer
        )
        cursor.execute(f"""
            UPDATE nodes n
            SET display_name = s.display_name,
                xrefs = s.xrefs
            FROM {staging_table} s
            WHERE n.node_id = s.node_id
        """)
        updated = cursor.rowcount
        cursor.execute(f"DROP TABLE IF EXISTS {staging_table}")
    raw_conn.commit()
    logger.info(f"Updated {updated} node(s) from {input_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--export', metavar='CSV_PATH', help="Write current node_id,display_name,xrefs to CSV_PATH.")
    group.add_argument('--apply', metavar='CSV_PATH', help="Read node_id,display_name,xrefs from CSV_PATH and UPDATE matching nodes in place.")
    parser.add_argument('--data-type', default=None, help="Only used with --export: restrict to nodes.data_type = this value (e.g. 'protein').")
    args = parser.parse_args()

    engine = _get_engine()
    SessionFactory = sessionmaker(bind=engine)
    session = SessionFactory()

    if args.export:
        export_nodes(session, args.export, args.data_type)
    else:
        apply_fixes(session, args.apply)

    session.close()
