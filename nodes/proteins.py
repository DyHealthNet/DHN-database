import json
import pandas as pd
import requests
from utils.query_nedrex import get_edge_associations
from nedrex.core import iter_nodes
from utils.models import Protein, ProteinAssocProtein
from settings import DEBUG


def get_protein_nodes(uniprot_ids: set[str] = None, observation_source: str = None) -> tuple[list[Protein], set[str]]:
    print("UniProt IDs:", len(uniprot_ids))
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
    print("Number of nans in UniProt col:", df['UniProt'].isna().sum())
    df['UniProt'] = df['UniProt'].fillna('')
    uniprot_ids = df['UniProt'].str.split('|')
    uniprot_ids = set([uniprot for sublist in uniprot_ids for uniprot in sublist])
    return uniprot_ids


def retrieve_interacting_proteins_neo4j(protein_ids):
    # Constructing a string of protein IDs for the Cypher query
    protein_id_string = ', '.join([f'"{protein_id}"' for protein_id in protein_ids])
    #P04049
    query = f"""
    MATCH (p1:Protein)-[r]->(p2:Protein)
    WHERE p1.primaryDomainId IN [{protein_id_string}] AND p2.primaryDomainId IN [{protein_id_string}]
    RETURN p1.primaryDomainId, p2.primaryDomainId, r
    """
    url = "https://api.nedrex.net/neo4j/query"
    response = requests.get(url, params={"query": query}, stream=True)
    #print("response")
    for line in response.iter_lines():
        print("response")
        print(json.loads(line.decode()))


def get_protein_interactions(proteinIds):
    # prefixed_proteinIds = {f"uniprot.{entry}" for entry in proteinIds}
    # retrieve_interacting_proteins_neo4j(proteinIds)
    assoc_graph = get_edge_associations(proteinIds, edge_type='protein_interacts_with_protein',
                                        direction='undirected')
    proteinInteractions = []
    for edge in assoc_graph.edges():
        uniprot_id_memberOne = edge[0]
        uniprot_id_memberTwo = edge[1]
        if (uniprot_id_memberTwo in proteinIds and uniprot_id_memberOne in proteinIds):
            proteinInteractions.append(ProteinAssocProtein(uniprot_id_1=uniprot_id_memberOne,
                                                           uniprot_id_2=uniprot_id_memberTwo))
    return proteinInteractions
