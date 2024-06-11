import os.path
import pandas as pd
import xml.etree.ElementTree as ET
from query_nedrex import domain_id_to_mondo, get_disorder_data


def download_metabolite_data(data_path: str):
    file_path = os.path.join(data_path, 'hmdb_metabolites.xml')
    if os.path.exists(file_path):
        return
    raise FileNotFoundError(f"Metabolite data at '{file_path}' does not exist, please download it from HMDB: "
                            f"https://hmdb.ca/downloads (All Metabolites in XML Format) and save the extracted file "
                            f"as 'hmdb_metabolites.xml' in the data directory")


def read_metabolite_mapping(mapping_file) -> pd.DataFrame:
    """
    Read the metabolite mapping file and return a dataframe
    """
    return pd.read_csv(mapping_file, sep='\t')


def read_hmdb_data(hmdb_file: str, relevant_ids: set[str], ext_ref: list = None) -> dict[str, dict]:
    """
    Read the HMDB data and return a dictionary of relevant metabolites
    :param hmdb_file: str, path to the HMDB file
    :param relevant_ids: set, set of relevant metabolite ids
    :param ext_ref: list, list of external references to extract, must end with '_id'
    """
    if ext_ref is None:
        ext_ref = ['kegg_id', 'chemspider_id', 'drugbank_id', 'pdb_id', 'wikipedia_id']

    assert all([True for ref in ext_ref if ref.endswith('_id')])

    tag = 'metabolite'
    context = ET.iterparse(hmdb_file, events=("start", "end"))
    context = iter(context)
    event, root = next(context)  # Get the root element
    hmdb_info = {}

    for event, elem in context:
        # remove the namespace
        elem.tag = elem.tag.split('}')[1] if '}' in elem.tag else elem.tag
        if not event == "end" or not elem.tag == tag:
            continue

        # get the accession as well as the secondary accession
        accession = elem.find('accession').text
        secondary_accessions = [sec.text for sec in elem.findall('secondary_accessions/accession')]
        if not accession in relevant_ids and not any(sec in relevant_ids for sec in secondary_accessions):
            root.clear()
            continue

        # get the display name, description, xrefs, synonyms and associated proteins
        display_name = elem.find('name').text
        description = elem.find('description').text
        # get xrefs from kegg, chemspider, drugbank, pdb, wikipedia
        xrefs = []
        for ref in ext_ref:
            xref = elem.find(f'{ref}')
            if xref is not None and xref.text is not None:
                xrefs.append(f"{ref.split('_')[0]}.{xref.text}")

        # add secondary accessions to xrefs
        xrefs.extend(secondary_accessions)
        synonyms = [synonym.text for synonym in elem.findall('synonyms/synonym')]
        proteins = [protein.text for protein in elem.findall('protein_association/protein/uniprot_id')]
        disease_associations = [disease.text for disease in elem.findall('diseases/disease/omim_id')
                                if disease.text is not None]
        hmdb_info[accession] = {'display_name': display_name,
                                'description': description,
                                'xrefs': xrefs,
                                'synonyms': synonyms,
                                'proteins': proteins,
                                'diseases': disease_associations}
        # clear the element
        root.clear()

    return hmdb_info


if __name__ == '__main__':
    metabo_data_path = '../data/DyHealthNet/chris_summary_data/metabolites/CHRIS_biocristes7500SumStats.txt'
    hmdb_data_path = '../data/hmdb_metabolites.xml'

