import os
import json
import pandas as pd
import networkx as nx
import urllib.request

def read_proteinID_chris(proteinID_path: str) :
    """
    reads Protein IDs from Chris dataset
    """
    df = pd.read_csv(proteinID_path, sep='\t')
    print(f"Found {len(df)} proteinIDs with {len(df['protein_id'].unique())} unique ids")
    return df

read_proteinID_chris("../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt")
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
