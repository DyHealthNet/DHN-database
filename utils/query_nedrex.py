import requests
import networkx as nx
import nedrex
from nedrex.core import api_keys_active, get_api_key
from utils.settings import DEBUG, ID_PREFIX
from nedrex.core import iter_nodes, iter_edges

# nedrex.config.set_url_base("https://api.nedrex.net/licensed/")
nedrex.config.set_url_base("https://apps.cosy.bio/licensed")
if api_keys_active():
    api_key = get_api_key(accept_eula=True)
    nedrex.config.set_api_key(api_key)


def get_disorder_data(pheno_ids: set[str]) -> list[dict]:
    """
    Fetches disorder data from nedrex for a set of phenotype reference ids
    :param pheno_ids: set of phenotype reference ids to fetch data for
    :return: list of dictionaries with disorder data
    """
    return [node for node in iter_nodes('disorder') if any(domain_id in pheno_ids for domain_id in node['domainIds'])]


def get_phenotype_data(hpo_ids: set[str]) -> list[dict]:
    """
    Fetches phenotype data from nedrex for a set of phenotype reference ids
    :param hpo_ids: set of phenotype reference ids to fetch data for
    :return: list of dictionaries with phenotype data
    """
    return [node for node in iter_nodes('phenotype') if node['primaryDomainId'] in hpo_ids]


def get_gene_data(entrez_ids: set[str] = None, search_col: str = 'primaryDomainId') -> list[dict]:
    """
    Fetches gene data from nedrex for a set of entrez ids
    :param entrez_ids: set of entrez ids to fetch data for
    :param search_col: column to search for entrez ids
    :return: list of dictionaries with gene data
    """
    if not entrez_ids:
        return [node for node in iter_nodes('gene')]
    return [node for node in iter_nodes('gene') if node[search_col] in entrez_ids]


def get_harmonizome_data(mondo_id: str) -> dict | None:
    """
    Fetches entrez ids associated with a mondo id as well as respective sources
    :param mondo_id: mondo id to fetch data for
    :return: dictionary with entrez ids and sources
    """
    # Deprecated API call
    return None
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
def domain_id_to_mondo(disorder_data: list, domain_id: str = ID_PREFIX) -> dict[str, str]:
    """
    Creates a dictionary with domain_id codes as keys and mondo ids as values
    :param disorder_data: dictionary with disorder data from nedrex
    :param domain_id: domain id to use for the mapping, i.e. snomedct, omim, orpha
    :return: dictionary with domain_id codes as keys and mondo ids as values
    """
    pheno_id_to_mondo = {}
    # go through all drug data and check if it has a domain_id code
    for disorder in disorder_data:
        if not 'domainIds' in disorder:
            continue

        for domain in disorder['domainIds']:
            if domain.startswith(domain_id):
                pheno_id_to_mondo[domain] = disorder['primaryDomainId']
    return pheno_id_to_mondo

#TODO why directed?
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

    if DEBUG:
        edges = []
        for edge in iter_edges(edge_type):
            if edge[first_node] in node_ids or edge[second_node] in node_ids:
                edges.append(edge)
            if len(edges) > 100:
                break
    else:
        edges = [e for e in iter_edges(edge_type) if e[first_node] in node_ids or e[second_node] in node_ids]

    G = nx.Graph()
    for association in edges:
        G.add_edge(association[first_node], association[second_node], source=association['dataSources'])
    return G
