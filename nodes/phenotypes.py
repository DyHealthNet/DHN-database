import os
import json
import pandas as pd
import networkx as nx
import urllib.request

import requests
from sqlalchemy import text
from sqlalchemy.orm import Session

from utils.query_nedrex import domain_id_to_mondo, get_edge_associations, get_harmonizome_data
from utils.settings import COHORT_COLUMNS, HPO_ID_PREFIX, ID_PREFIX, INPUT_ID_DB
from utils.models import Gene, Disorder, GeneAssocDisorder, Phenotype, DisorderAssocPhenotype
from utils.logger import get_logger

logger = get_logger(__name__)


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


def terms_from_hpo():
    """
    Retrieve all HPO terms from the HPO API
    :return: json data with all HPO terms
    """
    logger.debug("Downloading HPO terms from the HPO API")
    url = 'https://ontology.jax.org/api/hp/terms'
    response = requests.get(url)
    return response.json()


def read_hpo_ontology(hpo_path: str) -> dict:
    """
    Reads the HPO ontology json file
    :param hpo_path: path to the HPO ontology json file
    :return: dictionary with the HPO ontology data
    """
    with open(hpo_path, 'r') as f:
        hpo = json.load(f)
    return hpo


def get_needed_pheno_ids(phenotype_path: str) -> set[str]:
    """
    Extracts phenotype reference ids from a phenotype file
    :param phenotype_path: path to phenotype file
    :return: dataframe with phenotype data
    """
    file_name, ending = os.path.splitext(phenotype_path)
    if ending not in ['.csv', '.tsv']:
        raise ValueError(f"Unsupported file format for phenotypes meta file: {ending}. "
                         f"Only CSV and TSV files are supported.")
    sep = "," if ending == ".csv" else "\t"
    df = pd.read_csv(phenotype_path, sep=sep)
    unique_pheno = df[COHORT_COLUMNS['phenotype']['xref']].unique()
    pheno_ids = [pheno for sublist in [str(pheno).split(';') for pheno in unique_pheno] for pheno in sublist]
    logger.debug(f"Found {len(df)} phenotypes with {len(pheno_ids)} unique phenotype ids")
    pheno_ids = set([pheno.strip() for pheno in pheno_ids if pheno != "nan"])
    return pheno_ids


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


def pheno_ids_from_hpo_api(hpo_data: list, needed_pheno_ids: set[str]) -> dict:
    pheno_ids = {}
    for node in hpo_data:
        xrefs = node.get('xrefs', [])
        for ref in xrefs:
            if isinstance(HPO_ID_PREFIX, str):
                if HPO_ID_PREFIX not in ref:
                    continue
            elif not any(prefix in ref for prefix in HPO_ID_PREFIX):
                continue
            pheno_id = ref.split(':')[-1]
            if pheno_id in needed_pheno_ids:
                if pheno_id in pheno_ids:
                    pheno_ids[pheno_id].append(node['id'])
                else:
                    pheno_ids[pheno_id] = [node['id']]
    return pheno_ids


def pheno_ids_from_hpo(hpo_data: nx.Graph | list, needed_pheno_ids: set[str]) -> dict:
    if isinstance(hpo_data, nx.Graph):
        pheno_ids = pheno_ids_from_hpo_graph(hpo_data, needed_pheno_ids)
    else:
        pheno_ids = pheno_ids_from_hpo_api(hpo_data, needed_pheno_ids)
    return pheno_ids


def pheno_ids_from_hpo_graph(hpo_graph, needed_pheno_ids) -> dict:
    """
    Find the phenotype reference ids in the HPO ontology that are needed
    :param hpo_graph: Graph of the HPO ontology
    :param needed_pheno_ids: set of phenotype reference ids that are needed
    :return: mapping from phenotype reference ids to HPO ids
    """
    pheno_ids = {}
    for node in hpo_graph.nodes(data=True):
        xrefs = node[1].get('xrefs', [])
        for xref in xrefs:
            #TODO check if change back (delete if else and decomment)
            if isinstance(HPO_ID_PREFIX, str):
                if not xref['val'].startswith(HPO_ID_PREFIX):
                    continue
            elif not any(xref['val'].startswith(prefix) for prefix in HPO_ID_PREFIX):
                continue
            else:
            #if xref['val'].startswith(HPO_ID_PREFIX):
                pheno_id = xref['val'].split(':')[-1]
                if pheno_id in needed_pheno_ids:
                    hpo_id = node[0].split('/')[-1].replace('_', ':')
                    pheno_ids[pheno_id] = hpo_id
    return pheno_ids


#TODO return directly if current pheno_id are already OMIM/ORPHA Ids
def hpo_to_xref(data_dir, pheno_id_mapping) -> dict:
    """
    Reads the HPOA file and maps the phenotype reference ids to OMIM or ORPHA ids
    :param data_dir: directory where the HPOA file is stored
    :param pheno_id_mapping: mapping from phenotype reference ids to HPO ids
    :return: mapping from phenotype reference ids to OMIM or ORPHA ids
    """
    # read the HPOA file
    hpoa = pd.read_csv(f'{data_dir}/phenotype.hpoa', sep='\t', comment='#', low_memory=False)
    # convert the hpoa to a dict with hpo_id as key, database_id as value
    hpoa_database = hpoa.set_index('hpo_id')['database_id'].to_dict()
    # convert the pheno_id_mapping to a dict with pheno_id as key, hpo_id as value
    curr_pheno_id_to_omim = {}
    for key, value in pheno_id_mapping.items():
        try:
            curr_pheno_id_to_omim[key] = hpoa_database[value]
        except KeyError:
            pass
    return curr_pheno_id_to_omim

def disorder_to_mondo(disease_data, pheno_id_to_db) -> dict:
    """
    Convert the phenotype reference ids to OMIM or ORPHA ids and then to Mondo ids
    :param disease_data: disease data from the NEDREx API
    :param pheno_id_to_db: mapping from phenotype reference ids to OMIM or ORPHA ids
    :return: mapping from phenotype reference ids to Mondo ids
    """
    omim_ids = domain_id_to_mondo(disease_data, 'omim')
    orpha_ids = domain_id_to_mondo(disease_data, 'orpha')
    pheno_id_to_mondo = {}
    for key, value in pheno_id_to_db.items():
        raw_id = value.split(':')[-1]
        if value.startswith('OMIM'):
            pheno_id_to_mondo[key] = omim_ids.get(raw_id, None)
        elif value.startswith('ORPHA'):
            pheno_id_to_mondo[key] = orpha_ids.get(raw_id, None)
    # remove the None values
    return {k: v for k, v in pheno_id_to_mondo.items() if v is not None}


def mondo_to_phenotype(pheno_data, assoc_graph, pheno_id_hpo_map) -> dict:
    """
    Convert the Mondo ids to phenotype ids
    :param pheno_data: phenotype data from the NEDREx API
    :param assoc_graph: graph with associations between disorders and phenotypes
    :return: mapping from Mondo ids to phenotype (hpo) ids
    """
    mondo_to_pheno = {}
    relevant_hpo_ids = set(pheno_id_hpo_map.values())
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


def retrieve_disorder_data(needed_pheno_ids: set[str], pheno_ids_to_mondo: dict[str, str], descriptions: dict[str, str],
                           xrefs: dict, display_names: dict, gene_info: dict, assoc_graph: nx.Graph,
                           obs_source: str = None) -> tuple[set, set, set, set]:
    """
    Queries the needed phenotype reference ids and retrieves the associated genes and disorders from NeDRex
    :param display_names: Display names for the mondo ids
    :param gene_info: Information about the genes needed for the database (display name, synonyms, etc.)
    :param xrefs: cross references for the mondo ids to other databases
    :param descriptions: descriptions for the mondo ids
    :param needed_pheno_ids: phenotype reference ids in the dataset
    :param pheno_ids_to_mondo: map from phenotype reference ids to mondo ids from nedrex
    :param assoc_graph: association graph from nedrex of mondo ids to genes
    :param obs_source: Describes the source of observations - e.g. CHRIS
    :return: list of genes to add, list of disorders to add, list of gene associations to add,
    number of phenotype reference ids found
    """
    gene_associations = set()
    genes_to_add = set()
    disorders = set()
    #found = 0
    found = set()
    for pheno in needed_pheno_ids:
        # check if ; in phenotype reference id and if so, do this for all ids
        pheno_ids = str(pheno).split(';') #TODO this should already been taken care of in get_needed_pheno_ids() and if not Ids are probably wrong because the prefix is added element before
        for pheno_id in pheno_ids:
            mondo_id = pheno_ids_to_mondo.get(pheno_id, None)
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
            #found += 1
            found.add(pheno_id.split('.')[-1])
        continue
    return genes_to_add, disorders, gene_associations, found


def retrieve_phenotype_data(available_ids: dict, additional_data: dict, obs_source: str = None) \
        -> tuple[set, set, set, int]:
    """
    Retrieve the phenotype data for the needed phenotype reference ids
    # HPO conversion: Pathway
    # HPO data (HPO_ID ----> INPUT_ID) - look for needed phenotype reference IDs
    # -> map to Mondo (INPUT_ID -- OMIM_ID/ORPHA_ID --> Mondo_ID)
    #

    :param obs_source: Describes the source of observations - e.g. CHRIS
    :param additional_data: dictionary with additional data for the hpo ids, must be a dictionary with hpo ids as keys
    :param available_ids: dictionary with phenotype ids as keys and hpo ids as values
    :return: dictionary with the phenotype data
    """
    found = 0
    genes_to_add = set() #TODO why is this here
    phenotypes = set()
    disorder_associations = set()

    available_pheno_ids = available_ids
    # go through all the nodes in the HPO graph and find the ones that have xrefs to the phenotype ID DB #TODO?

    logger.debug(f'Found {len(available_pheno_ids)} {INPUT_ID_DB} ids in the HPO ontology')

    # get edge associations for disorder_has_phenotype
    assoc_graph = get_edge_associations(set(el for v in available_pheno_ids.values() for el in v), edge_type='disorder_has_phenotype')

    # find the genes that are associated with the mondo ids
    for pheno_id, hpo_ids in available_pheno_ids.items():
        pheno_id = f"{ID_PREFIX}.{pheno_id}"
        for hpo_id in hpo_ids:
            if additional_data.get(hpo_id, None) is None:
                continue
            phenotype_data = additional_data[hpo_id]
            phenotypes.add(Phenotype(hpo_id=hpo_id,
                                     xrefs=set(phenotype_data['domainIds'] + [pheno_id]),
                                     description=phenotype_data['description'],
                                     synonyms=phenotype_data['synonyms'],
                                     display_name=phenotype_data['displayName'],
                                     observation_source=obs_source))
            found += 1
            if hpo_id not in assoc_graph:
                continue
            # get the disorder ids associated with the hpo id
            for edge in assoc_graph.edges(hpo_id, data=True):
                disorder = edge[1]
                source = edge[2]['source'][0]
                new_assoc = DisorderAssocPhenotype(mondo_id=disorder, hpo_id=hpo_id, edge_source=source)
                disorder_associations.add(new_assoc)

    return genes_to_add, phenotypes, disorder_associations, found


def get_additional_diseases(session: Session, obs_source: str = None):
    # I know this defeats the purpose of SQLAlchemy but I could not find a way to do this with the ORM
    sql_string = f"""SELECT xrefs
                    FROM disorder
                    WHERE EXISTS (
                        SELECT 1
                        FROM unnest(xrefs) AS xref
                        WHERE xref LIKE 'omim.%'
                    ) AND observation_source = '{obs_source}';"""
    return {x for x in session.execute(text(sql_string)).fetchall() for x in x[0] if x.startswith('omim.')}
