import os
import json
import pandas as pd
import networkx as nx
import urllib.request

import requests
from sqlalchemy import text

from utils.query_nedrex import domain_id_to_mondo, get_edge_associations, \
    get_harmonizome_data, get_phenotype_data

from utils.models import Gene, Disorder, GeneAssocDisorder, Phenotype, DisorderAssocPhenotype


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


def get_needed_snomed_ids(phenotype_path: str) -> set[str]:
    """
    Extracts snomed ids from a phenotype file
    :param phenotype_path: path to phenotype file
    :return: dataframe with phenotype data
    """
    df = pd.read_csv(phenotype_path, sep='\t')
    uniqe_snomed = df['snomed_id'].unique()
    snomeds = [snomed for sublist in [str(snomed).split(';') for snomed in uniqe_snomed] for snomed in sublist]
    print(f"Found {len(df)} phenotypes with {len(snomeds)} unique snomed ids")
    snomeds = set([f"snomedct.{snomed}" for snomed in snomeds if snomed != "nan"])
    return snomeds


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
    needed_snomed_ids = set([x.split(".")[1] for x in needed_snomed_ids])
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


def pheno_pathway(available_snomed_ids, assoc_graph):
    # go the phenotype way through nedrex
    # first get the available hpo ids from nedrex
    # phenotype_data = is mapping from hpo ids to mondo ids
    phenotype_assoc_graph = get_edge_associations(edge='disorder_has_phenotype')
    available_snomed_ids = {snomed_id: hpo_id.replace(':', '.').replace('HP', 'hpo') for snomed_id, hpo_id in
                            available_snomed_ids.items()}

    phenotype_data = domain_id_to_mondo(get_phenotype_data(available_snomed_ids), 'hpo')

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


def mondo_in_association_graph(mondo_id: str, assoc_graph: nx.Graph) -> tuple[list[str], list[str]] | tuple[None, None]:
    """
    Retrieves the genes associated with a mondo id from the association graph
    :param mondo_id: A mondo id, can be None
    :param assoc_graph: The association graph from NEDRex
    :return: genes associated with the mondo id, source of the data, harmonizome data if not in the association graph
    """
    if mondo_id is None:
        return None, None
    if mondo_id in assoc_graph:
        genes = []
        sources = []
        for edge in assoc_graph.edges(mondo_id, data=True):
            genes.append(edge[1])
            sources.append(edge[2]['source'][0])
    else:
        harm = get_harmonizome_data(mondo_id)
        if harm is None:
            return None, None
        genes = []
        sources = []
        for key, value in harm.items():
            genes.extend(value)
            sources.extend([key] * len(value))
    return genes, sources


def retrieve_disorder_data(needed_snomed: set[str], snomed_to_mondo: dict[str, str], descriptions: dict[str, str],
                           xrefs: dict, display_names: dict, gene_info: dict, assoc_graph: nx.Graph,
                           obs_source: str = None) -> tuple[set, set, set, int]:
    """
    Queries the needed snomed ids and retrieves the associated genes and disorders from NeDRex
    :param display_names: Display names for the mondo ids
    :param gene_info: Information about the genes needed for the database (display name, synonyms, etc.)
    :param xrefs: cross references for the mondo ids to other databases
    :param descriptions: descriptions for the mondo ids
    :param needed_snomed: snomed ids in the dataset
    :param snomed_to_mondo: map from snomed to mondo ids from nedrex
    :param assoc_graph: association graph from nedrex of mondo ids to genes
    :param obs_source: Describes the source of observations - e.g. CHRIS
    :return: list of genes to add, list of disorders to add, list of gene associations to add,
    number of snomed ids found
    """
    gene_associations = set()
    genes_to_add = set()
    disorders = set()
    found = 0
    for snomed in needed_snomed:
        # check if ; in snomed id and if so, do this for all ids
        snomed_ids = str(snomed).split(';')
        for snomed_id in snomed_ids:
            mondo_id = snomed_to_mondo.get(snomed_id)
            xref = xrefs.get(mondo_id, None)
            description = descriptions.get(mondo_id, None)
            display_name = display_names.get(mondo_id, None)
            genes, sources = mondo_in_association_graph(mondo_id, assoc_graph)
            if genes is None:
                continue
            # add genes to set
            for gene in genes:
                gene_data = gene_info.get(gene, None)
                if gene_data is None:
                    new_gene = Gene(entrez_id=gene, observation_source='external')
                    genes_to_add.add(new_gene)
                    continue
                new_gene = Gene(entrez_id=gene_data['primaryDomainId'],
                                display_name=gene_data['displayName'],
                                description=gene_data['description'],
                                synonyms=gene_data['synonyms'],
                                chromosome=gene_data['chromosome'],
                                observation_source='external')
                genes_to_add.add(new_gene)

            disorders.add(
                Disorder(mondo_id=mondo_id,
                         display_name=display_name,
                         xrefs=xref,
                         description=description,
                         observation_source=obs_source))
            # add gene associations to set for each source
            gene_associations.update([GeneAssocDisorder(entrez_id=gene, mondo_id=mondo_id, edge_source=source)
                                      for gene, source in zip(genes, sources)])
            found += 1
        continue
    return genes_to_add, disorders, gene_associations, found


def retrieve_phenotype_data(available_ids: dict, additional_data: dict, obs_source: str = None) \
        -> tuple[set, set, set, int]:
    """
    Retrieve the phenotype data for the needed snomed ids
    # HPO conversion: Pathway
    # HPO data (HPO_ID ----> SNOMED_ID) - look for needed SNOMED IDs
    # -> map to Mondo (SNOMED_ID -- OMIM_ID/ORPHA_ID --> Mondo_ID)
    #

    :param obs_source: Describes the source of observations - e.g. CHRIS
    :param additional_data: dictionary with additional data for the hpo ids, must be a dictionary with hpo ids as keys
    :param available_ids: dictionary with snomed ids as keys and hpo ids as values
    :return: dictionary with the phenotype data
    """
    found = 0
    genes_to_add = set()
    phenotypes = set()
    disorder_associations = set()

    available_snomed_ids = available_ids
    # go through all the nodes in the HPO graph and find the ones that have xrefs to SNOMED

    print(f'Found {len(available_snomed_ids)} snomed ids in the HPO ontology')

    # get edge associations for disorder_has_phenotype
    assoc_graph = get_edge_associations(set(available_snomed_ids.values()), edge_type='disorder_has_phenotype')

    # find the genes that are associated with the mondo ids
    for snomed_id, hpo_id in available_snomed_ids.items():
        snomed_id = f"snomedct.{snomed_id}"
        if additional_data.get(hpo_id, None) is None:
            continue
        phenotype_data = additional_data[hpo_id]
        phenotypes.add(Phenotype(hpo_id=hpo_id,
                                 xrefs=set(phenotype_data['domainIds'] + [snomed_id]),
                                 description=phenotype_data['description'],
                                 synonyms=phenotype_data['synonyms'],
                                 display_name=phenotype_data['displayName'],
                                 observation_source=obs_source))
        if hpo_id not in assoc_graph:
            continue
        # get the disorder ids associated with the hpo id
        for edge in assoc_graph.edges(hpo_id, data=True):
            disorder = edge[1]
            source = edge[2]['source'][0]
            new_assoc = DisorderAssocPhenotype(mondo_id=disorder, hpo_id=hpo_id, edge_source=source)
            disorder_associations.add(new_assoc)

        found += 1
    return genes_to_add, phenotypes, disorder_associations, found


def get_additional_diseases(session, obs_source: str = None):
    # I know this defeats the purpose of SQLAlchemy but I could not find a way to do this with the ORM
    sql_string = f"""SELECT xrefs
                    FROM disorder
                    WHERE EXISTS (
                        SELECT 1
                        FROM unnest(xrefs) AS xref
                        WHERE xref LIKE 'omim.%'
                    ) AND observation_source = '{obs_source}';"""
    return {x for x in session.execute(text(sql_string)).fetchall() for x in x[0] if x.startswith('omim.')}


if __name__ == '__main__':

    # data handling
    data_dir = '../../data'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    needed_files = [f'{data_dir}/hp.json', f'{data_dir}/phenotype.hpoa']
    if not all([os.path.exists(f) for f in needed_files]):
        download_hpo_ontology(data_dir)

    hpo_data = read_hpo_ontology(needed_files[0])
    hpo_graph = ontology_data_to_network(hpo_data)
    needed_ids = get_needed_snomed_ids('../../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv')

    # go through all the nodes in the HPO graph and find the ones that have xrefs to SNOMED
    available_snomed_ids = snomed_from_hpo(hpo_graph, needed_ids)
