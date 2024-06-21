# load the calculated edges from the file and add them to the database
# steps:
# 1. load the calculated edges from the file
# 2. load info files for label names
# 3. convert label names to proper names
# 4. add the edges to the database

import pandas as pd
from models import *

DB_EDGES = {
    ('protein', 'protein'): (EffectsProteinProtein, 'uniprot_id_1', 'uniprot_id_2'),
    ('protein', 'phenotype'): (EffectsProteinPhenotype, 'uniprot_id', 'hpo_id'),
    ('protein', 'disorder'): (EffectsProteinMetabolite, 'uniport_id', 'mondo_id'),
    ('protein', 'metabolite'): (EffectsProteinMetabolite, 'uniprot_id', 'hmdb_id'),

    ('metabolite', 'metabolite'): (EffectsMetaboliteMetabolite, 'hmdb_id_1', 'hmdb_id_2'),
    ('metabolite', 'disorder'): (EffectsMetaboliteDisorder, 'hmdb_id', 'mondo_id'),
    ('metabolite', 'phenotype'): (EffectsMetabolitePhenotype, 'hmdb_id', 'hpo_id'),

    ('phenotype', 'phenotype'): (EffectsPhenotypePhenotype, 'hpo_id_1', 'hpo_id_2'),
    ('phenotype', 'disorder'): (EffectsPhenotypeDisorder, 'hpo_id', 'mondo_id'),

    ('disorder', 'disorder'): (EffectsDisorderDisorder, 'mondo_id_1', 'mondo_id_2')

}


def load_files(file_path, sep="\t"):
    return pd.read_csv(file_path, sep=sep)


def get_labels(label_file: pd.DataFrame, label_type: str = None):
    name_col = {
        'protein': ('protein_id', 'UniProt'),
        'phenotype': ('label', 'snomed_id'),
        'metabolite': (None, 'hmdb_id')
    }
    label_col, id_col = name_col[label_type]
    # check that the label and id cols exist
    assert label_col in label_file.columns, f"Label column {label_col} not found in the file"
    assert id_col in label_file.columns, f"ID column {id_col} not found in the file"

    return label_file.set_index(label_col)[id_col].to_dict()


def diff_phenotype_disorder():
    pass


def map_edge(edge: tuple[str, str], protein_map: dict, pheno_map: dict, metabo_map: dict):
    source_type = None
    target_type = None
    mapped_source = None
    mapped_target = None
    for map, type in zip([protein_map, pheno_map, metabo_map], ['protein', 'phenotype', 'metabolite']):
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


def format_edges(edges: pd.DataFrame, protein_map: dict, phenotype_map: dict, metabolite_map: dict):
    formatted_edges = []
    for edge_row in edges.iterrows():
        edge_row = edge_row[1]
        edge = edge_row['label1'], edge_row['label2']
        mapped_edge, types = map_edge(edge, protein_map, phenotype_map, metabolite_map)
        if not types:
            print(f"Edge {edge} not found in any of the maps")
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
    return formatted_edges


if __name__ == '__main__':
    edges_path = '../data/scores.csv'
    pheno_data_path = '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv'
    protein_data_path = '../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt'
    metabo_data_path = '../data/DyHealthNet/chris_summary_data/metabolites/CHRIS_biocristes7500SumStats.txt'

    phenotypes = load_files(pheno_data_path)
    proteins = load_files(protein_data_path)
    metabolites = load_files(metabo_data_path)
    edges = load_files(edges_path, sep=",")

    pheno_map = get_labels(phenotypes, 'phenotype')
    protein_map = get_labels(proteins, 'protein')
    try:
        metabo_map = get_labels(metabolites, 'metabolite')
    except AssertionError:
        metabo_map = None

    formatted_edges = format_edges(edges, protein_map, pheno_map, metabo_map)
