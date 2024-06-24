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


DB_EDGES = {
    ('protein', 'protein'): (EffectsProteinProtein, 'uniprot_id_1', 'uniprot_id_2'),
    ('protein', 'phenotype'): (EffectsProteinPhenotype, 'uniprot_id', 'hpo_id'),
    ('phenotype', 'protein'): (EffectsProteinPhenotype, 'hpo_id', 'uniprot_id'),
    ('protein', 'disorder'): (EffectsProteinDisorder, 'uniprot_id', 'mondo_id'),
    ('disorder', 'protein'): (EffectsProteinDisorder, 'mondo_id', 'uniprot_id'),
    ('protein', 'metabolite'): (EffectsProteinMetabolite, 'uniprot_id', 'hmdb_id'),
    ('metabolite', 'protein'): (EffectsProteinMetabolite, 'hmdb_id', 'uniprot_id'),

    ('metabolite', 'metabolite'): (EffectsMetaboliteMetabolite, 'hmdb_id_1', 'hmdb_id_2'),
    ('metabolite', 'disorder'): (EffectsMetaboliteDisorder, 'hmdb_id', 'mondo_id'),
    ('disorder', 'metabolite'): (EffectsMetaboliteDisorder, 'mondo_id', 'hmdb_id'),
    ('metabolite', 'phenotype'): (EffectsMetabolitePhenotype, 'hmdb_id', 'hpo_id'),
    ('phenotype', 'metabolite'): (EffectsMetabolitePhenotype, 'hpo_id', 'hmdb_id'),

    ('phenotype', 'phenotype'): (EffectsPhenotypePhenotype, 'hpo_id_1', 'hpo_id_2'),
    ('phenotype', 'disorder'): (EffectsPhenotypeDisorder, 'hpo_id', 'mondo_id'),
    ('disorder', 'phenotype'): (EffectsPhenotypeDisorder, 'mondo_id', 'hpo_id'),

    ('disorder', 'disorder'): (EffectsDisorderDisorder, 'mondo_id_1', 'mondo_id_2')

}


def load_files(file_path, sep="\t"):
    return pd.read_csv(file_path, sep=sep)


def get_labels(label_file: pd.DataFrame, label_type: str = None):
    name_col = {
        'protein': ('protein_id', 'UniProt'),
        'phenotype': ('label', 'snomed_id'),
        'metabolite': ('analyte_name', 'hmdb_id')
    }
    label_col, id_col = name_col[label_type]
    # check that the label and id cols exist
    assert label_col in label_file.columns, f"Label column {label_col} not found in the file"
    assert id_col in label_file.columns, f"ID column {id_col} not found in the file"

    return label_file.set_index(label_col)[id_col].to_dict()


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


def map_edge(edge: tuple[str, str], protein_map: dict, pheno_map: dict, metabo_map: dict, disorder_map: dict = None):
    source_type = None
    target_type = None
    mapped_source = None
    mapped_target = None
    for map, type in zip([protein_map, pheno_map, metabo_map, disorder_map],
                         ['protein', 'phenotype', 'metabolite', 'disorder']):
        # map the source and target to the proper names
        if not mapped_source:
            mapped_source = map.get(edge[0], None)
            if mapped_source:
                source_type = type

        if not mapped_target:
            mapped_target = map.get(edge[1], None)
            if mapped_target:
                target_type = type

        if source_type and target_type:
            break

    if not source_type or not target_type:
        return None, None
    return (mapped_source, mapped_target), (source_type, target_type)


def format_edges(edges: pd.DataFrame, protein_map: dict, phenotype_map: dict, metabolite_map: dict, disorder_map: dict):
    formatted_edges = []
    num_edge_types = {}
    for idx, edge_row in tqdm.tqdm(edges.iterrows(), maxinterval=len(edges), desc="Rows"):
        edge = edge_row['label1'], edge_row['label2']
        mapped_edge, types = map_edge(edge, protein_map, phenotype_map, metabolite_map, disorder_map)
        if not types:
            # print(f"Edge {edge} not found in any of the maps")
            continue
        source, target = mapped_edge
        source_type, target_type = types
        edge_type, source_col, target_col = DB_EDGES[(source_type, target_type)]
        edge_values = {source_col: source,
                       target_col: target,
                       'p_value': edge_row['pval'],
                       'adjusted_p_value': edge_row['adj_pval'],
                       'effect_size': edge_row['effsize'],
                       'effect_size_type': edge_row['effsize_type']}
        edge = edge_type(**edge_values)
        formatted_edges.append(edge)
        num_edge_types[edge_type] = num_edge_types.get(edge_type, 0) + 1
    for edge_type, count in num_edge_types.items():
        print(f"Added {count} edges of type {edge_type}")
    return formatted_edges


def add_edges(session, edges):
    session.add_all(edges)
    try:
        session.commit()
    except Exception as e:
        session.rollback()
        print(f"A problem occurred while adding edges: {e}")
        return False


def main(session, edges_path, pheno_data_path, protein_data_path, metabo_data_path):
    phenotypes = load_files(pheno_data_path)
    proteins = load_files(protein_data_path)
    metabolites = load_files(metabo_data_path)
    edges = load_files(edges_path, sep=",")

    # get base labels
    pheno_base_map = get_labels(phenotypes, 'phenotype')
    protein_map = get_labels(proteins, 'protein')
    metabo_map = get_labels(metabolites, 'metabolite')

    # format the labels and filter out the ones that do not exist in the database
    existing_proteins = set(filter_exising_ids(session, Protein.uniprot_id))
    existing_metabolites = set(filter_exising_ids(session, Metabolite.hmdb_id))
    protein_map = {k: f"uniprot.{v}" for k, v in protein_map.items() if v in existing_proteins}
    metabo_map = {k: f"hmdb.{v}" for k, v in metabo_map.items() if v in existing_metabolites}
    pheno_map, disorder_map = diff_phenotype_disorder(pheno_base_map, session)

    # formatted_edges = format_edges(edges, protein_map, pheno_map, metabo_map, disorder_map)
    formatted_edges = format_edges(edges[edges['pval'] <= 0.05], protein_map, pheno_map, metabo_map, disorder_map)
    add_edges(session, formatted_edges)


if __name__ == '__main__':
    url = url_object = URL.create(
        "postgresql",
        username="postgres",
        password="password",  # plain (unescaped) text
        host="0.0.0.0",
        port=9000,
        database="postgres",
    )
    engine = create_engine(url)
    Session = sessionmaker(bind=engine)
    session = Session()

    edges_path = '../data/scores.csv'
    pheno_data_path = '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv'
    protein_data_path = '../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt'
    metabo_data_path = '../data/DyHealthNet/chris_summary_data/metabolites/CHRIS_biocristes7500SumStats.txt'

    main(session, edges_path, pheno_data_path, protein_data_path, metabo_data_path)

