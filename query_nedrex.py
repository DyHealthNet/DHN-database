import requests


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
def snomedct_to_mondo(disorder_data: dict) -> dict:
    """
    Creates a dictionary with snomedct codes as keys and mondo ids as values
    :param disorder_data: dictionary with disorder data from nedrex
    :return: dictionary with snomedct codes as keys and mondo ids as values
    """
    snomed_to_mondo = {}
    # go through all drug data and check if it has a snomedct code
    for disorder in disorder_data:
        if not 'domainIds' in disorder:
            continue

        for domain in disorder['domainIds']:
            if domain.startswith('snomedct'):
                snomed_id = domain.split('.')[1]
                snomed_to_mondo[snomed_id] = disorder['primaryDomainId']
    return snomed_to_mondo


if __name__ == '__main__':
    data = get_disorder_data()
    snomed_to_mondo = snomedct_to_mondo(data)
    # pick random snomedct code
    snomed_id = list(snomed_to_mondo.keys())[0]
    print(f"Snomedct code: {snomed_id}")
    harmonizome_data = get_harmonizome_data(snomed_to_mondo[snomed_id])
    print(harmonizome_data)
