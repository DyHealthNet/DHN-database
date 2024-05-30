import os
import json
import pandas as pd
import networkx as nx
import urllib.request

import requests

from query_nedrex import get_disorder_data, domain_id_to_mondo, needed_snomed_ids, get_edge_associations, \
    get_harmonizome_data, get_phenotype_data


def download_hpo_ontology(data_dir: str) -> None:
    """
    Download the HPO ontology files from the latest release on github
    :param data_dir: directory where the files should be stored
    :return: None
    """
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


def read_hpo_ontology(hpo_path: str) -> dict:
    """
    Reads the HPO ontology json file
    :param hpo_path: path to the HPO ontology json file
    :return: dictionary with the HPO ontology data
    """
    with open(hpo_path, 'r') as f:
        hpo = json.load(f)
    return hpo


def ontology_data_to_network(hpo_data: dict) -> nx.Graph:
    """
    Convert the HPO ontology data to a networkx graph
    :param hpo_data: downloaded HPO ontology data
    :return: graph of the HPO ontology with nodes as HPO ids and edges as relationships
    """
    graph = nx.Graph()
    for node in hpo_data['graphs'][0]['nodes']:
        node_id = node['id']
        xrefs = node.get('meta', {}).get('xrefs', [])
        graph.add_node(node_id, xrefs=xrefs)

    for edge in hpo_data['graphs'][0]['edges']:
        graph.add_edge(edge['sub'], edge['obj'], edge=edge['pred'])
    return graph


def snomed_from_hpo(hpo_graph, needed_snomed_ids) -> dict:
    """
    Find the SNOMED ids in the HPO ontology that are needed
    :param hpo_graph: Graph of the HPO ontology
    :param needed_snomed_ids: set of snomed ids that are needed
    :return: mapping from snomed ids to HPO ids
    """
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


def hpo_to_xref(data_dir, snomed_id_mapping) -> dict:
    """
    Reads the HPOA file and maps the snomed ids to OMIM or ORPHA ids
    :param data_dir: directory where the HPOA file is stored
    :param snomed_id_mapping: mapping from snomed ids to HPO ids
    :return: mapping from snomed ids to OMIM or ORPHA ids
    """
    # read the HPOA file
    hpoa = pd.read_csv(f'{data_dir}/phenotype.hpoa', sep='\t', comment='#', low_memory=False)
    # convert the hpoa to a dict with hpo_id as key, database_id as value
    hpoa_database = hpoa.set_index('hpo_id')['database_id'].to_dict()
    # convert the snomed_id_mapping to a dict with snomed_id as key, hpo_id as value
    snomed_to_omim = {}
    for key, value in snomed_id_mapping.items():
        try:
            snomed_to_omim[key] = hpoa_database[value]
        except KeyError:
            pass
    return snomed_to_omim


def disorder_to_mondo(disease_data, snomed_to_db) -> dict:
    """
    Convert the snomed ids to OMIM or ORPHA ids and then to Mondo ids
    :param disease_data: disease data from the NEDREx API
    :param snomed_to_db: mapping from snomed ids to OMIM or ORPHA ids
    :return: mapping from snomed ids to Mondo ids
    """
    omim_ids = domain_id_to_mondo(disease_data, 'omim')
    orpha_ids = domain_id_to_mondo(disease_data, 'orpha')
    snomed_to_mondo = {}
    for key, value in snomed_to_db.items():
        raw_id = value.split(':')[-1]
        if value.startswith('OMIM'):
            snomed_to_mondo[key] = omim_ids.get(raw_id, None)
        elif value.startswith('ORPHA'):
            snomed_to_mondo[key] = orpha_ids.get(raw_id, None)
    # remove the None values
    return {k: v for k, v in snomed_to_mondo.items() if v is not None}


def mondo_to_phenotype(pheno_data, assoc_graph, snomed_hpo_map) -> dict:
    """
    Convert the Mondo ids to phenotype ids
    :param pheno_data: phenotype data from the NEDREx API
    :param assoc_graph: graph with associations between disorders and phenotypes
    :return: mapping from Mondo ids to phenotype (hpo) ids
    """
    mondo_to_pheno = {}
    relevant_hpo_ids = set(snomed_hpo_map.values())
    for _, pheno_id in pheno_data.items():
        if pheno_id in assoc_graph and pheno_id in relevant_hpo_ids:
            mondo_ids = assoc_graph[pheno_id]
            if len(mondo_ids) == 1:
                try:
                    mondo_to_pheno[mondo_ids[0]] = [pheno_id]
                except KeyError:
                    continue
                continue

            for m_id in mondo_ids:
                try:
                    mondo_to_pheno[m_id].append(pheno_id)
                except KeyError:
                    mondo_to_pheno[m_id] = [pheno_id]
    return mondo_to_pheno


def omim_pathway(available_snomed_ids, data_dir, assoc_graph):
    # convert the snomed ids to OMIM ids
    snomed_to_xref = hpo_to_xref(data_dir, available_snomed_ids)
    print(f'Found {len(snomed_to_xref)} snomed ids with OMIM/ORPHA ids')

    # get the disorder data
    disorder_data = get_disorder_data()
    final_mapping = disorder_to_mondo(disorder_data, snomed_to_xref)
    print(f'Found {len(final_mapping)} snomed ids with Mondo ids')
    # find the genes that are associated with the mondo ids
    found = 0
    found_ids = {}
    for snomed_id, mondo_id in final_mapping.items():
        if mondo_id in assoc_graph or get_harmonizome_data(mondo_id):
            found += 1
            found_ids[snomed_id] = mondo_id
    print(f'Found {found} snomed ids with associated genes')
    return found_ids


def pheno_pathway(available_snomed_ids, assoc_graph):
    # go the phenotype way through nedrex
    # first get the available hpo ids from nedrex
    phenotype_data = domain_id_to_mondo(get_phenotype_data(), 'hpo')
    # get the associations between phenotypes and disorders (mondo ids - hpo ids)
    phenotype_assoc_graph = get_edge_associations(edge='disorder_has_phenotype')
    available_snomed_ids = {snomed_id: hpo_id.replace(':', '.').replace('HP', 'hpo') for snomed_id, hpo_id in
                            available_snomed_ids.items()}

    # map the mondo ids to the phenotype ids
    phenotype_mapping = mondo_to_phenotype(phenotype_data, phenotype_assoc_graph, available_snomed_ids)

    # now for every mondo id check its associated genes and map them to the hpo + snomed ids
    found_pheno_ids = {}
    for mondo_id, pheno_id in phenotype_mapping.items():
        if mondo_id in assoc_graph or get_harmonizome_data(mondo_id):
            # convert every pheno_id to snomed ids
            snomed_ids = [snomed for hpo_id in pheno_id for snomed, hpo in available_snomed_ids.items() if hpo == hpo_id]
            for snomed_id in snomed_ids:
                found_pheno_ids[snomed_id] = mondo_id
    print(f'Found {len(found_pheno_ids)} snomed ids with associated genes through phenotypes')
    return found_pheno_ids


if __name__ == '__main__':

    # data handling
    data_dir = '../data'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    needed_files = [f'{data_dir}/hp.json', f'{data_dir}/phenotype.hpoa']
    if not all([os.path.exists(f) for f in needed_files]):
        download_hpo_ontology(data_dir)

    hpo_data = read_hpo_ontology(needed_files[0])
    hpo_graph = ontology_data_to_network(hpo_data)
    needed_ids = needed_snomed_ids('../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv')
    needed_ids = set(needed_ids['snomed_id'].unique())
    assoc_graph = get_edge_associations()

    # go through all the nodes in the HPO graph and find the ones that have xrefs to SNOMED
    available_snomed_ids = snomed_from_hpo(hpo_graph, needed_ids)

    print(f'Found {len(available_snomed_ids)} snomed ids in the HPO ontology')
    found_ids = set(omim_pathway(available_snomed_ids, data_dir, assoc_graph))
    found_pheno_ids = set(pheno_pathway(available_snomed_ids, assoc_graph))
    print(f'Ovelapping snomed ids: {len(set(found_ids).intersection(found_pheno_ids))}')