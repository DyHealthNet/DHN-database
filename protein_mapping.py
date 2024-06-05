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



def get_edge_associations2(node_ids: set[str], edge_type) -> nx.Graph:
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
    edges = iter_edges(edge_type)
    #edges = [e for e in iter_edges(edge_type) if e[sourceDomainId] in node_ids or e[targetDomainId] in node_ids]
    G = nx.Graph()
    for association in edges:
        G.add_edge(association[sourceDomainId], association[targetDomainId], source=association['dataSources'])
    return G



def add_proteinSet_data(session, protein_data_path):
    proteinData = read_proteinID_chris(protein_data_path)
    proteinSet = set(proteinData)
    associationsType = list["protein_encoded_by",
                            "protein_has_signature",
                            "protein_in_pathway",
                            "protein_interacts_with_protein",
                            "protein_similarity_protein"]
    protein_encoded_by_graph = get_edge_associations2(node_ids=proteinSet, edge_type='protein_interacts_with_protein')
    for edge in protein_encoded_by_graph:
        protein = edge[1]
        source = edge[2]['source'][0] #frage hier ist das eigentlichj das Objekt Protein mit den Attributen und ich kann hier sequence fetchen?
        newProt = Protein()

    #disorder_associations.add(new_assoc)




#proteinIDs_CHRIS = read_proteinID_chris("../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt")

#for proteinID in proteinIDs_CHRIS:
 #   #print(proteinID)
  #  print(get_proteinID_neddrex(proteinID)[0]['geneName'])

["is_isoform_of","molecule_similarity_molecule","protein_encoded_by","protein_has_signature","protein_in_pathway","protein_interacts_with_protein","protein_similarity_protein"]




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
