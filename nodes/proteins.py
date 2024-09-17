import json
import pandas as pd
import requests
from utils.query_nedrex import get_edge_associations
from nedrex.core import iter_nodes
from utils.models import Protein, ProteinAssocProtein
from utils.logger import get_logger
from settings import DEBUG

logger = get_logger(__name__)


def get_protein_nodes(uniprot_ids: set[str] = None, observation_source: str = None) -> tuple[list[Protein], set[str]]:
    logger.debug(f"UniProt IDs: {len(uniprot_ids)}")
    protein_set = []
    found_proteins = set()
    uniprot_ids = {f"uniprot.{uniprot_id}" for uniprot_id in uniprot_ids}
    for node in iter_nodes('protein'):
        # Remove the 'uniprot.' prefix from node primaryDomainId
        primary_domain_id = node['primaryDomainId']
        if primary_domain_id in uniprot_ids:
            protein = Protein(
                uniprot_id=str(primary_domain_id),
                display_name=f"{node.get('displayName').split('_')[0]}",  # Remove the species from the display name
                gene_entrez_id=str(node.get('geneName')),
                sequence=str(node.get('sequence')),
                description=str(node.get('comments')),
                observation_source=observation_source
            )
            protein_set.append(protein)
            found_proteins.add(primary_domain_id)
        if DEBUG:
            if len(protein_set) > 100:
                break

    return protein_set, found_proteins


def read_protein_id_chris(proteinID_path: str):
    """
    reads Protein IDs from Chris dataset
    """
    df = pd.read_csv(proteinID_path, sep='\t')
    # check how many nans in the uniprot ocl
    logger.debug(f"Number of nans in UniProt col: {df['UniProt'].isna().sum()}")
    df['UniProt'] = df['UniProt'].fillna('')
    uniprot_ids = df['UniProt'].str.split('|')
    uniprot_ids = set([uniprot for sublist in uniprot_ids for uniprot in sublist])
    return uniprot_ids


def get_protein_interactions(proteinIds):
    # prefixed_proteinIds = {f"uniprot.{entry}" for entry in proteinIds}
    # retrieve_interacting_proteins_neo4j(proteinIds)
    assoc_graph = get_edge_associations(proteinIds, edge_type='protein_interacts_with_protein',
                                        direction='undirected')
    protein_interactions = []
    for edge in assoc_graph.edges():
        uniprot_id_member_one = edge[0]
        uniprot_id_member_two = edge[1]
        if uniprot_id_member_two in proteinIds and uniprot_id_member_one in proteinIds:
            protein_interactions.append(ProteinAssocProtein(uniprot_id_1=uniprot_id_member_one,
                                                            uniprot_id_2=uniprot_id_member_two))
    return protein_interactions
