import os
import json
import pandas as pd
import networkx as nx
import urllib.request

import requests

from query_nedrex import get_disorder_data, domain_id_to_mondo, needed_snomed_ids


def download_hpo_ontology(data_dir: str):
    # find the latest release of the HPO ontology
    url = 'https://api.github.com/repos/obophenotype/human-phenotype-ontology/releases/latest'
    response = requests.get(url)
    response = response.json()
    # check for the assets in the release
    mapping_link, phenotype_link = None, None
    for asset in response['assets']:
        if asset['name'] == 'hp.json':
            mapping_link = asset['browser_download_url']
        if asset['name'] == 'phenotype.hpoa':
            phenotype_link = asset['browser_download_url']
    if not mapping_link or not phenotype_link:
        raise ValueError('Could not find the HPO ontology files, please check the release on github and download '
                         'the hp.json and phenotype.hpoa files manually')

    download_path = f'{data_dir}/hp.json'
    urllib.request.urlretrieve(mapping_link, download_path)
    download_path = f'{data_dir}/phenotype.hpoa'
    urllib.request.urlretrieve(phenotype_link, download_path)


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


def snomed_from_hpo(hpo_graph, needed_snomed_ids):
    snomed_ids = {}
    for node in hpo_graph.nodes(data=True):
        xrefs = node[1].get('xrefs', [])
        for xref in xrefs:
            if xref['val'].startswith('SNOMEDCT_US:'):
                snomed_id = xref['val'].split(':')[-1]
                if snomed_id in needed_snomed_ids:
                    hpo_id = node[0].split('/')[-1].replace('_', ':')
                    snomed_ids[snomed_id] = hpo_id
    return snomed_ids


def hpo_to_omim(data_dir, snomed_id_mapping):
    # read the HPOA file
    hpoa = pd.read_csv(f'{data_dir}/phenotype.hpoa', sep='\t', comment='#', low_memory=False)
    # convert the hpoa to a dict with hpo_id as key, database_id as value
    hpoa_dict = hpoa.set_index('hpo_id')['database_id'].to_dict()
    # convert the snomed_id_mapping to a dict with snomed_id as key, hpo_id as value
    snomed_to_omim = {}
    for key, value in snomed_id_mapping.items():
        try:
            snomed_to_omim[key] = hpoa_dict[value]
        except KeyError:
            pass
    return snomed_to_omim


if __name__ == '__main__':
    # data handling
    data_dir = '../data'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    needed_files = [f'{data_dir}/hp.json', f'{data_dir}/phenotype.hpoa']
    if not all([os.path.exists(f) for f in needed_files]):
        download_hpo_ontology(data_dir)

    # HPO conversion
    hpo_data = read_hpo_ontology(needed_files[0])
    hpo_graph = ontology_data_to_network(hpo_data)
    needed_ids = needed_snomed_ids('../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv')
    needed_ids = set(needed_ids['snomed_id'].unique())
    # go through all the nodes in the HPO graph and find the ones that have xrefs to SNOMED
    available_snomed_ids = snomed_from_hpo(hpo_graph, needed_ids)

    # convert the snomed ids to OMIM ids
    snomed_to_omim = hpo_to_omim(data_dir, available_snomed_ids)
    print(f'Found {len(snomed_to_omim)} snomed ids with OMIM ids')
    print(f'Found {len(available_snomed_ids)} snomed ids in the HPO ontology')



