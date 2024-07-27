# load the calculated edges from the file and add them to the database
# steps:
# 1. load the calculated edges from the file
# 2. load info files for label names
# 3. convert label names to proper names
# 4. add the edges to the database
import timeit

import pandas as pd
import tqdm

from models import *
from sqlalchemy import URL, create_engine, text
from sqlalchemy.orm import sessionmaker
from settings import DEBUG

DB_EDGES = {
    ('protein', 'protein'): (EffectsProteinProtein, 'protein_id_1', 'protein_id_2'),
    ('protein', 'phenotype'): (EffectsProteinPhenotype, 'protein_id', 'phenotype_id'),
    ('phenotype', 'protein'): (EffectsProteinPhenotype, 'phenotype_id', 'protein_id'),
    ('protein', 'metabolite'): (EffectsProteinMetabolite, 'protein_id', 'metabolite_id'),
    ('metabolite', 'protein'): (EffectsProteinMetabolite, 'metabolite_id', 'protein_id'),

    ('metabolite', 'metabolite'): (EffectsMetaboliteMetabolite, 'metabolite_id_1', 'metabolite_id_2'),
    ('metabolite', 'phenotype'): (EffectsMetabolitePhenotype, 'metabolite_id', 'phenotype_id'),
    ('phenotype', 'metabolite'): (EffectsMetabolitePhenotype, 'phenotype_id', 'metabolite_id'),

    ('phenotype', 'phenotype'): (EffectsPhenotypePhenotype, 'phenotype_id_1', 'phenotype_id_2'),
}


def load_files(file_path, sep="\t"):
    return pd.read_csv(file_path, sep=sep)


def get_labels(label_file: pd.DataFrame, label_type: str = None):
    name_col = {
        'protein': ('protein_id', 'protein_id'),
        'metabolite': ('analyte_name', 'analyte_name'),
        'phenotype': ('label', 'label')
    }
    label_col, id_col = name_col[label_type]
    # check that the label and id cols exist
    assert label_col in label_file.columns, f"Label column {label_col} not found in the file"
    assert id_col in label_file.columns, f"ID column {id_col} not found in the file"

    identity_dict = label_file[label_col].to_dict()
    return {value: value for key, value in identity_dict.items()}


def get_xref_rows(session, table: str = 'disorders', xref_col: str = 'omim'):
    # I know this defeats the purpose of SQLAlchemy but I could not find a way to do this with the ORM
    sql_string = f"""SELECT *
                    FROM {table}
                    WHERE EXISTS (
                        SELECT 1
                        FROM unnest(xrefs) AS xref
                        WHERE xref LIKE '{xref_col}.%'
                    );"""
    return session.execute(text(sql_string)).fetchall()


def filter_exising_ids(session, column):
    return session.query(column).all()


def diff_phenotype_disorder(base_map: dict, session):
    phenotypes = get_xref_rows(session, 'phenotypes', 'snomedct')
    snomed_hpo_map = {[xref.split(".")[1] for xref in x[3] if xref.startswith('snomedct')][0]: x[0] for x in phenotypes}
    disorders = get_xref_rows(session, 'disorders', 'snomedct')
    snomed_mondo_map = {[xref.split(".")[1] for xref in x[2] if xref.startswith('snomedct')][0]: x[0] for x in disorders}
    phenotype_map = dict()
    disorder_map = dict()
    for key, val in base_map.items():
        if val in snomed_hpo_map:
            phenotype_map[key] = snomed_hpo_map[val]
        elif val in snomed_mondo_map:
            disorder_map[key] = snomed_mondo_map[val]
    return phenotype_map, disorder_map


def map_edge(edge: pd.Series, protein_map: dict, pheno_map: dict, metabo_map: dict):
    maps_and_types = [(protein_map, 'protein'), (pheno_map, 'phenotype'), (metabo_map, 'metabolite')]

    mapped_source = mapped_target = source_type = target_type = None

    for map, type in maps_and_types:
        if not mapped_source:
            mapped_source = map.get(edge['label1'])
            if mapped_source:
                source_type = type

        if not mapped_target:
            mapped_target = map.get(edge['label2'])
            if mapped_target:
                target_type = type

        if source_type and target_type:
            break

    if not source_type or not target_type:
        return None, None
    return (mapped_source, mapped_target), (source_type, target_type)


def process_chunk(edges_chunk, protein_map, phenotype_map, metabolite_map):

    def map_and_filter(edge):
        mapped, types = map_edge(edge, protein_map, phenotype_map, metabolite_map)
        return mapped, types if types else None

    # Apply the map_and_filter function to the chunk
    mapped_results = edges_chunk.apply(map_and_filter, axis=1)

    # Filter the valid results
    valid_results = [result for result in mapped_results if result[1] is not None]

    if not valid_results:
        return [], []

    # Separate the mapped edges and types
    mapped_edges, types_list = zip(*valid_results)

    # Extract valid rows based on indices of valid results
    valid_indices = [i for i, result in enumerate(mapped_results) if result[1] is not None]
    valid_edges = edges_chunk.iloc[valid_indices]

    # Create a new DataFrame with the mapped data
    formatted_edges = pd.DataFrame({
        'source': [source for source, _ in mapped_edges],
        'target': [target for _, target in mapped_edges],
        'source_type': [source_type for source_type, _ in types_list],
        'target_type': [target_type for _, target_type in types_list],
        'p_value': valid_edges['pval'].values,
        'adjusted_p_value': valid_edges['adj_pval'].values,
        'effect_size': valid_edges['effsize'].values,
        'effect_size_type': valid_edges['effsize_type'].values
    })

    def create_edge_data(row):
        edge_type, source_col, target_col = DB_EDGES[(row['source_type'], row['target_type'])]
        return {
            'edge_type': edge_type,
            source_col: row['source'],
            target_col: row['target'],
            'p_value': row['p_value'],
            'adjusted_p_value': row['adjusted_p_value'],
            'effect_size': row['effect_size'],
            'effect_size_type': row['effect_size_type']
        }

    # Apply the create_edge_data function
    formatted_edges['edge_data'] = formatted_edges.apply(create_edge_data, axis=1)

    # Create the final formatted edges list
    formatted_edges_list = [
        row['edge_type'](**{key: value for key, value in row.items() if key != 'edge_type'})
        for row in formatted_edges['edge_data']
    ]

    return formatted_edges_list, [edge.__class__ for edge in formatted_edges_list]


def format_edges(session, edges: pd.DataFrame, protein_map: dict, phenotype_map: dict, metabolite_map: dict):
    all_edge_types = {}
    num_edge_types = {}
    chunk_size = 10_000_000
    num_chunks = (len(edges) // chunk_size) + 1

    if DEBUG:
        chunk_size = 10_000
        num_chunks = 1
        edges = edges.sample(frac=1)

    for i in range(num_chunks):
        start_index = i * chunk_size
        end_index = min((i + 1) * chunk_size, len(edges))
        edges_chunk = edges.iloc[start_index:end_index]

        # Process the current chunk
        formatted_edges_list, edge_types = process_chunk(edges_chunk, protein_map, phenotype_map, metabolite_map)

        add_edges(session, formatted_edges_list)
        del formatted_edges_list
        print(f"Chunk {i + 1}/{num_chunks} added successfully")
        # do the value counts of the edges and add them to a running total
        chunk_edge_types = pd.Series(edge_types).value_counts().to_dict()
        num_edge_types = {edge_type: num_edge_types + chunk_edge_types.get(edge_type, 0)
                          for edge_type, num_edge_types in all_edge_types.items()}

    for edge_type, count in num_edge_types.items():
        print(f"Added {count} edges of type {edge_type}")
    return


def add_edges(session, edges):
    session.add_all(edges)
    try:
        session.commit()
    except Exception as e:
        session.rollback()
        print(f"A problem occurred while adding edges: {e}")
        return False


def add_calculated_edges(session, edges_path, pheno_data_path, protein_data_path, metabo_data_path):
    phenotypes = load_files(pheno_data_path)
    proteins = load_files(protein_data_path)
    metabolites = load_files(metabo_data_path)
    edges = load_files(edges_path, sep=",")

    # get base labels
    pheno_map = get_labels(phenotypes, 'phenotype')
    protein_map = get_labels(proteins, 'protein')
    metabo_map = get_labels(metabolites, 'metabolite')
    print("All maps loaded successfully")

    # pheno_map, disorder_map = diff_phenotype_disorder(pheno_base_map, session)

    format_edges(session, edges, protein_map, pheno_map, metabo_map)
    # formatted_edges = format_edges(edges[edges['pval'] <= 0.05], protein_map, pheno_map, metabo_map, disorder_map)


if __name__ == '__main__':
    url = url_object = URL.create(
        "postgresql",
        username="postgres",
        password="password",  # plain (unescaped) text
        host="0.0.0.0",
        port=9852,
        database="postgres",
    )
    engine = create_engine(url)
    Session = sessionmaker(bind=engine)
    session = Session()

    edges_path = '../data/scores.csv'
    pheno_data_path = '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv'
    protein_data_path = '../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt'
    metabo_data_path = '../data/DyHealthNet/chris_summary_data/metabolites/CHRIS_biocristes7500SumStats.txt'
    add_calculated_edges(session, edges_path, pheno_data_path, protein_data_path, metabo_data_path)
