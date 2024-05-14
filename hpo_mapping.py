import os
import json
import networkx as nx
import urllib.request
from query_nedrex import get_disorder_data, domain_id_to_mondo, needed_snomed_ids


def download_hpo_ontology(data_dir: str):
    # find the latest release of the HPO ontology
    download_link = 'https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2024-04-26/hp.json'
    download_path = f'{data_dir}/hp.json'
    urllib.request.urlretrieve(download_link, download_path)


def read_hpo_ontology(hpo_path: str):
    with open(hpo_path, 'r') as f:
        hpo = json.load(f)
    return hpo


def ontology_data_to_network(hpo_data: dict):
    graph = nx.Graph()
    for node in hpo_data['graphs'][0]['nodes']:
        node_id = node['id']
        xrefs = node.get('meta', {}).get('xrefs', [])
        graph.add_node(node_id, xrefs=xrefs)

    for edge in hpo_data['graphs'][0]['edges']:
        graph.add_edge(edge['sub'], edge['obj'], edge=edge['pred'])
    return graph


if __name__ == '__main__':
    data_dir = '../data'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    hpo_path = f'{data_dir}/hp.json'
    if not os.path.exists(hpo_path):
        download_hpo_ontology(data_dir)
    hpo_data = read_hpo_ontology(hpo_path)
    hpo_graph = ontology_data_to_network(hpo_data)
    needed_ids = needed_snomed_ids('../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv')
    needed_ids = set(needed_ids['snomed_id'].unique())
    # go through all the nodes in the HPO graph and find the ones that have xrefs to SNOMED
    available_snomed_ids = set()
    for node in hpo_graph.nodes(data=True):
        xrefs = node[1].get('xrefs', [])
        for xref in xrefs:
            if xref['val'].startswith('SNOMEDCT_US:'):
                snomed_id = xref['val'].split(':')[-1]
                available_snomed_ids.add(snomed_id)
    print(f'Found {len(available_snomed_ids)} snomed ids in the HPO ontology')
    print(f'Found {len(needed_ids)} snomed ids in the DyHealthNet data')
    print(f'Found {len(needed_ids.intersection(available_snomed_ids))} snomed ids in both')

