import requests
import networkx as nx
import pandas as pd


def get_disorder_data() -> dict:
    """
    Fetches all disorder data from nedrex
    :return: dictionary with disorder data
    """
    url = 'https://api.nedrex.net/disorder/all'
    response = requests.get(url)
    data = response.json()
    return data


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
    except:
        return None
    return data


# this function should be in another file
def domain_id_to_mondo(disorder_data: dict, domain_id: str = 'snomedct') -> dict:
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
                snomed_id = domain.split('.')[1]
                snomed_to_mondo[snomed_id] = disorder['primaryDomainId']
    return snomed_to_mondo


def get_all_associations() -> nx.Graph:
    """
    Fetches all associations from NedRex
    :return: networkx graph with all mondo ids and associated genes
    """
    url = 'https://api.nedrex.net/open/gene_associated_with_disorder/all'
    response = requests.get(url)
    data = response.json()
    G = nx.Graph()
    types = set()
    for association in data:
        G.add_edge(association['sourceDomainId'], association['targetDomainId'])
        types.add(association['type'])
    return G


# this function should be in another file
def needed_snomed_ids(phenotype_path: str):
    """
    Extracts snomed ids from a phenotype file
    :param phenotype_path: path to phenotype file
    :return: dataframe with phenotype data
    """
    df = pd.read_csv(phenotype_path, sep='\t')
    print(f"Found {len(df)} phenotypes with {len(df['snomed_id'].unique())} unique snomed ids")
    return df


if __name__ == '__main__':
    needed_snomed_ids = needed_snomed_ids('../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv')
    data = get_disorder_data()
    snomed_to_mondo = domain_id_to_mondo(data)
    assoc_graph = get_all_associations()
    found = 0
    # go through all snomed ids and check if they have a mondo id
    for snomed_id in needed_snomed_ids['snomed_id'].unique():
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
