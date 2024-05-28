import os
import json
import pandas as pd
import networkx as nx
import urllib.request
import requests
import itertools



def read_proteinID_chris(proteinID_path: str) :
    """
    reads Protein IDs from Chris dataset
    """
    df = pd.read_csv(proteinID_path, sep='\t')
   # print(f"Found {len(df)} proteinIDs with {len(df['UniProt'].unique())} unique ids")
    #return df['UniProt'].unique()

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




#proteinIDs_CHRIS = read_proteinID_chris("../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt")

#for proteinID in proteinIDs_CHRIS:
 #   #print(proteinID)
  #  print(get_proteinID_neddrex(proteinID)[0]['geneName'])





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
