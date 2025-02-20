import os

import pandas as pd
from utils.query_nedrex import get_edge_associations, get_gene_data
from nedrex.core import iter_nodes
from utils.models import Protein, ProteinAssocProtein, ProteinAssocGene, Gene
from utils.logger import get_logger
from utils.settings import DEBUG

logger = get_logger(__name__)


def get_protein_nodes(uniprot_ids: set[str] = None, observation_source: str = None) \
        -> tuple[set[Protein], set[tuple]]:
    logger.debug(f"UniProt IDs: {len(uniprot_ids)}")
    proteins = set()
    protein_assoc_genes = set()
    uniprot_ids = {f"uniprot.{uniprot_id}" for uniprot_id in uniprot_ids}
    for node in iter_nodes('protein'):
        # Remove the 'uniprot.' prefix from node primaryDomainId
        primary_domain_id = node['primaryDomainId']
        if primary_domain_id in uniprot_ids:
            protein = Protein(
                uniprot_id=primary_domain_id,
                display_name=f"{node.get('displayName').split('_')[0]}",  # Remove the species from the display name
                sequence=node.get('sequence'),
                description=node.get('comments'),
                observation_source=observation_source
            )
            protein_assoc_genes.add((primary_domain_id, node.get('geneName')))
            proteins.add(protein)
        if DEBUG:
            if len(proteins) > 100:
                break

    return proteins, protein_assoc_genes


def read_protein_id_chris(protein_id_path: str):
    """
    reads Protein IDs from Chris dataset
    """
    file_name, ending = os.path.splitext(protein_id_path)
    if ending not in ['.csv', '.tsv']:
        raise ValueError(f"Unsupported file format for phenotypes meta file: {ending}. "
                         f"Only CSV and TSV files are supported.")
    sep = "," if ending == ".csv" else "\t"
    df = pd.read_csv(protein_id_path, sep=sep)
    # check how many nans in the uniprot ocl
    logger.debug(f"Number of nans in UniProt col: {df['UniProt'].isna().sum()}")
    df['UniProt'] = df['UniProt'].fillna('')
    uniprot_ids = df['UniProt'].str.split('|')
    uniprot_ids = set([uniprot for sublist in uniprot_ids for uniprot in sublist])
    return uniprot_ids


def get_protein_interactions(protein_ids):
    # prefixed_proteinIds = {f"uniprot.{entry}" for entry in proteinIds}
    # retrieve_interacting_proteins_neo4j(proteinIds)
    assoc_graph = get_edge_associations(protein_ids, edge_type='protein_interacts_with_protein',
                                        direction='undirected')
    protein_interactions = []
    for edge in assoc_graph.edges():
        uniprot_id_member_one = edge[0]
        uniprot_id_member_two = edge[1]
        if uniprot_id_member_two in protein_ids and uniprot_id_member_one in protein_ids:
            protein_interactions.append(ProteinAssocProtein(uniprot_id_1=uniprot_id_member_one,
                                                            uniprot_id_2=uniprot_id_member_two))
    return protein_interactions
