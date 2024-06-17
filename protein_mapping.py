import os
import json
import pandas as pd
import networkx as nx
import urllib.request
import requests
import itertools
from query_nedrex import get_disorder_data, domain_id_to_mondo, get_needed_snomed_ids, get_edge_associations, \
    get_harmonizome_data, get_phenotype_data
from nedrex.core import iter_nodes, iter_edges
from nedrex.core import api_keys_active, get_api_key




def read_proteinID_chris(proteinID_path: str) :
    """
    reads Protein IDs from Chris dataset
    """
    df = pd.read_csv(proteinID_path, sep='\t')
   # print(f"Found {len(df)} proteinIDs with {len(df['UniProt'].unique())} unique ids")

    return df['UniProt'][1:5].unique()

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



def get_edge_associations2(node_ids: set[str], edge_type):
    """
    Fetches all edges of a certain type that are associated with a set of node ids
    :param edge_type: type of edge to fetch
    :param node_ids: set of node ids to fetch edges for (i.e. mondo ids)
    :return: networkx graph with all mondo ids and associated genes
    """
    sourceDomainId = "sourceDomainId"
    targetDomainId = "targetDomainId"
    if(edge_type == "protein_interacts_with_protein"):
        sourceDomainId = "memberOne"
        targetDomainId = "memberTwo"

        edges = [e for e in iter_edges(edge_type) if e['memberOne'] in node_ids or e['memberTwo'] in node_ids]
        G = nx.Graph()
        for association in edges:
            G.add_edge(association['memberOne'], association['memberTwo'], source=association['dataSources'])
        return G



    #edgesGenerator = iter_edges(edge_type)
    #testdict = dict(edgesGenerator)
    #test =2
    #nodeGenerator = iter_nodes("Proteins")
   # test = [for node_ids in iter_edges]
    #filter(edges)
    #edges = [e for e in iter_edges(edge_type) if e[sourceDomainId] in node_ids or e[targetDomainId] in node_ids]
    #G = nx.Graph()
   # for association in edges:
    #    G.add_edge(association[sourceDomainId], association[targetDomainId], source=association['dataSources'])
    filtered_edges = []
    for edge in edgesGenerator:
        # Check if either the source or target node ID of the edge is in the node_ids set
        if edge[sourceDomainId] in node_ids or edge[targetDomainId] in node_ids:
            # If so, add the edge to the filtered list
            filtered_edges.append(edge)
    return edgesGenerator



def add_proteinSet_data(session, protein_data_path):
    proteinData = read_proteinID_chris(protein_data_path)
    proteinSet = set(proteinData)
    associationsType = list["protein_encoded_by",
                            "protein_has_signature",
                            "protein_in_pathway",
                            "protein_interacts_with_protein",
                            "protein_similarity_protein"]


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



import json
import requests
def testneo4j():
    query = "MATCH (n) RETURN n LIMIT 25"
    url = "http://nedrex-api.zbh.uni-hamburg.de/neo4j/query"
    response = requests.get(url, params={"query":query}, stream=True)
    for line in response.iter_lines():
        print(json.loads(line.decode()))