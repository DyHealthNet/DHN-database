import requests
import networkx as nx
import pandas as pd
import nedrex
from nedrex.core import iter_nodes, iter_edges
from nedrex.core import api_keys_active, get_api_key


#nedrex.config.set_url_base("https://api.nedrex.net/open/")
nedrex.config.set_url_base(" https://apps.cosy.bio/licensed")
if api_keys_active():
    api_key = get_api_key(accept_eula=True)
    nedrex.config.set_api_key(api_key)


def get_disorder_data(snomedct_ids: set[str]) -> list[dict]:
    """
    Fetches disorder data from nedrex for a set of snomedct ids
    :param snomedct_ids: set of snomedct ids to fetch data for
    :return: list of dictionaries with disorder data
    """
    return [node for node in iter_nodes('disorder') if any(domain_id in snomedct_ids for domain_id in node['domainIds'])]


def get_phenotype_data(hpo_ids: set[str]) -> list[dict]:
    """
    Fetches phenotype data from nedrex for a set of snomedct ids
    :param hpo_ids: set of snomedct ids to fetch data for
    :return: list of dictionaries with phenotype data
    """
    return [node for node in iter_nodes('phenotype') if node['primaryDomainId'] in hpo_ids]


def get_gene_data(entrez_ids: set[str] = None) -> list[dict]:
    """
    Fetches gene data from nedrex for a set of entrez ids
    :param entrez_ids: set of entrez ids to fetch data for
    :return: list of dictionaries with gene data
    """
    if not entrez_ids:
        return [node for node in iter_nodes('gene')]
    return [node for node in iter_nodes('gene') if node['primaryDomainId'] in entrez_ids]


def get_harmonizome_data(mondo_id: str) -> dict | None:
    """
    Fetches entrez ids associated with a mondo id as well as respective sources
    :param mondo_id: mondo id to fetch data for
    :return: dictionary with entrez ids and sources
    """
    url = f'https://api.nedrex.net/static/harmonizome/{mondo_id}'
    response = requests.get(url)
    try:
        data = response.json()
        source_gene_map = {}
        for gene_object in data:
            source = gene_object['source']
            source_gene_map[source] = gene_object['genes']
    except TypeError:
        return None
    return source_gene_map


# this function should be in another file
def domain_id_to_mondo(disorder_data: list, domain_id: str = 'snomedct') -> dict[str, str]:
    """
    Creates a dictionary with snomedct codes as keys and mondo ids as values
    :param disorder_data: dictionary with disorder data from nedrex
    :param domain_id: domain id to use for the mapping, i.e. snomedct, omim, orpha
    :return: dictionary with snomedct codes as keys and mondo ids as values
    """
    snomed_to_mondo = {}
    # go through all drug data and check if it has a snomedct code
    for disorder in disorder_data:
        if not 'domainIds' in disorder:
            continue

        for domain in disorder['domainIds']:
            if domain.startswith(domain_id):
                snomed_to_mondo[domain] = disorder['primaryDomainId']
    return snomed_to_mondo


def get_edge_associations(node_ids: set[str], edge_type='gene_associated_with_disorder', direction='directed') -> nx.Graph:
    """
    Fetches all edges of a certain type that are associated with a set of node ids
    :param direction: The direction of the edges to fetch, either 'directed' or 'undirected'
    :param edge_type: type of edge to fetch
    :param node_ids: set of node ids to fetch edges for (i.e. mondo ids)
    :return: networkx graph with all mondo ids and associated genes
    """
    if direction == 'directed':
        first_node = 'sourceDomainId'
        second_node = 'targetDomainId'
    elif direction == 'undirected':
        first_node = 'memberOne'
        second_node = 'memberTwo'
    else:
        raise ValueError(f"Direction {direction} not supported")

    # edges = [e for e in iter_edges(edge_type) if e[first_node] in node_ids or e[second_node] in node_ids]

    edges =[]
    count = 0
    for edge in iter_edges(edge_type):
        if count >= 2:
            break
        if edge[first_node] in node_ids or edge[second_node] in node_ids:
            edges.append(edge)
            count += 1

    G = nx.Graph()
    for association in edges:
        G.add_edge(association[first_node], association[second_node], source=association['dataSources'])
    return G


# this function should be in another file
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


if __name__ == '__main__':
    needed_snomeds = get_needed_snomed_ids('../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv')
    data = get_disorder_data(needed_snomeds)
    snomed_to_mondo = domain_id_to_mondo(data)
    assoc_graph = get_edge_associations(set(snomed_to_mondo.values()))
    found = 0

    # go through all snomed ids and check if they have a mondo id
    for snomed_id in needed_snomeds:
        # check if ; in snomed id and if so, do this for all ids
        snomed_ids = str(snomed_id).split(';')
        for snomed_id in snomed_ids:
            mondo_id = snomed_to_mondo.get(snomed_id)
            if mondo_id is None:
                print(f"No mondo id found for snomedct code {snomed_id}")
                continue
            if not mondo_id in assoc_graph and not get_harmonizome_data(mondo_id):
                print(f"Neither mondo id nor harmonizome data found for mondo id {mondo_id}")
                continue
            found += 1
        continue
    print(f"Found and successfully matched {found} snomed ids with diseases")
