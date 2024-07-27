import os
import json
import requests
import pandas as pd
import networkx as nx
import urllib.request
import requests
import itertools
from query_nedrex import get_disorder_data, domain_id_to_mondo, get_needed_snomed_ids, get_edge_associations, \
    get_harmonizome_data, get_phenotype_data
from nedrex.core import iter_nodes, iter_edges
from models import Protein
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
                uniprot_id= str(primary_domain_id),
                gene_entrez_id=str(node.get('geneName')),
                sequence=str(node.get('sequence')),
                description=str(node.get('comments')),
                observation_source=observation_source
            )
            protein_set.append(protein)
            found_proteins.add(primary_domain_id)
        if DEBUG:
            if len(protein_set) > 10:
                break

    return protein_set, found_proteins


def get_genomic_variants(entrez_ids: set[str] = None, observation_source: str = None) -> list[dict]:
    print("len(ids):", len(entrez_ids))
    variant_affects_gene_graph = get_edge_associations(node_ids=entrez_ids_list, edge_type='variant_affects_gene',
                                                       direction='directed')  # Graph with 1511628 nodes and 1539719 edges
    test = 1
    # print(entrez_ids)
    # entrez_id = variant_affects_gene_graph['entrez_id']
    # variant_primaryDomainId = variant_affects_gene_graph['variant_primaryDomainId']
    print(type(variant_affects_gene_graph))
    variant_primaryDomainId_list = []
    entrez_id_list = []
    for edge in variant_affects_gene_graph.edges(data=True):
        variant_primaryDomainId = edge[0]
        entrez_id = edge[1]
    test = 2
    protein_set = []
    for node in iter_edges('variant_affects_gene'):
        # Remove the 'uniprot.' prefix from node primaryDomainId
        primary_domain_id = node['primaryDomainId']
        if primary_domain_id in uniprot_ids:
            protein = Protein(
                uniprot_id= str(primary_domain_id),
                gene_entrez_id=str(node.get('geneName')),
                sequence=str(node.get('sequence')),
                description=str(node.get('comments')),
                observation_source=observation_source
            )
            protein_set.append(protein)
    return protein_set



def read_proteinID_chris(proteinID_path: str) :
    """
    reads Protein IDs from Chris dataset
    """
    df = pd.read_csv(proteinID_path, sep='\t')
    # check how many nans in the uniprot ocl
    print("Number of nans in UniProt col:", df['UniProt'].isna().sum())
    return df['UniProt'].unique()


def retrieve_interacting_proteins_neo4j(protein_ids):
    # Constructing a string of protein IDs for the Cypher query
    protein_id_string = ', '.join([f'"{protein_id}"' for protein_id in protein_ids])
    #P04049
    query= f"""
    MATCH (p1:Protein)-[r]->(p2:Protein)
    WHERE p1.primaryDomainId IN [{protein_id_string}] AND p2.primaryDomainId IN [{protein_id_string}]
    RETURN p1.primaryDomainId, p2.primaryDomainId, r
    """
    url = "https://api.nedrex.net/neo4j/query"
    response = requests.get(url, params={"query":query}, stream=True)
    #print("response")
    for line in response.iter_lines():
        print("response")
        print(json.loads(line.decode()))


#["is_isoform_of","molecule_similarity_molecule","protein_encoded_by","protein_has_signature","protein_in_pathway","protein_interacts_with_protein","protein_similarity_protein"]

"""
"protein_encoded_by",
"protein_has_signature",
"protein_in_pathway",
"protein_interacts_with_protein",
"protein_similarity_protein"
"""

[
  "disorder",
  "drug",
  "gene",
  "pathway",
  "protein",
  "signature"
]




def testneo4j():
    query = "MATCH (n) RETURN n LIMIT 25"
    url = "http://nedrex-api.zbh.uni-hamburg.de/neo4j/query"
    response = requests.get(url, params={"query":query}, stream=True)
    for line in response.iter_lines():
        print(json.loads(line.decode()))
def get_proteinID_neddrex(proteinID: str):
    """
    fetches data for proteinID from neddrex
    https://api.nedrex.net/relations/get_drugs_targetting_proteins
    https://api.nedrex.net/relations/get_encoded_proteins
    https://api.nedrex.net/get_by_id/protein?q=uniprot.Q9UBT6
    :param proteinID:
    :return:
    """
    url = f'https://api.nedrex.net/get_by_id/protein?q=uniprot.' + proteinID
    response = requests.get(url)
    try:
        data = response.json()
    except:
        return None
    return data
