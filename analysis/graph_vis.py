from matplotlib import pyplot as plt
from sqlalchemy import create_engine, inspect, URL, text, MetaData, Table
from sqlalchemy.orm import sessionmaker
from sqlalchemy import inspect as sql_inspect
import networkx as nx
from settings import *


def generate_graph(engine):
    edge_tables = ['effects', 'assoc', 'references']
    layer_1_tables = 'effects'
    layer_2_tables = 'assoc'
    reference_tables = 'references'
    inspector = sql_inspect(engine)
    G = nx.Graph()
    node_colors = []
    edge_labels = {}

    metadata = MetaData()
    metadata.reflect(bind=engine)

    # Add nodes (tables)
    for table_name in inspector.get_table_names():
        table = Table(table_name, metadata, autoload_with=engine)
        row_count = session.query(table).count()
        if row_count == 0:
            continue

        # check if substring of edge tables contains edge substrings
        if any([x in table_name for x in edge_tables]):
            continue
        G.add_node(table_name)
        if 'cohort' in table_name:
            node_colors.append('lightblue')
        else:
            node_colors.append('lightgreen')

    # Add edges (foreign key relationships)
    for table_name in inspector.get_table_names():
        parts = table_name.split("_")
        if layer_1_tables in table_name:
            node_1 = f"cohort_{parts[1]}"
            node_2 = f"cohort_{parts[2]}"
        elif reference_tables in table_name:
            node_1 = f"cohort_{parts[2]}"
            node_2 = f"{parts[2]}"
        elif layer_2_tables in table_name:
            node_1 = f"{parts[0]}"
            node_2 = f"{parts[2]}"
        else:
            continue
        # handle special case:
        if table_name == 'cohort_references_disease':
            node_1 = 'cohort_phenotype'
            node_2 = 'disorder'
        elif table_name == 'cohort_references_variant':
            node_1 = 'cohort_variant'
            node_2 = 'genomic_variant'
        elif table_name == 'variant_associates_gene':
            node_1 = 'genomic_variant'
            node_2 = 'gene'
        if not G.has_node(node_1) or not G.has_node(node_2):
            continue
        G.add_edge(node_1, node_2)
        edge_labels[(node_1, node_2)] = table_name

    return G, node_colors, edge_labels


url = url_object = URL.create(
    "postgresql",
    username=DB_USER,
    password=DB_PASSWORD,
    host=DB_HOST,
    port=DB_PORT,
    database=DB_NAME,
)

engine = create_engine(url)

Session = sessionmaker(bind=engine)
session = Session()

G, node_colors, edge_labels = generate_graph(engine)

assert len(G.nodes) == len(node_colors)

pos = nx.spring_layout(G, 0.7)
plt.figure(figsize=(10, 8), dpi=300)
nx.draw(G, pos, with_labels=True, node_size=3000, node_color=node_colors, font_size=10, font_weight='bold', edge_color='gray')
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='blue', font_size=8)
plt.savefig("uml_diagram.png")
plt.show()
