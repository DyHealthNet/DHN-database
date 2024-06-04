import pandas as pd


def read_metabolite_mapping(mapping_file):
    """
    Read the metabolite mapping file and return a dataframe
    """
    return pd.read_csv(mapping_file, sep='\t')


def query_hmdb():
    """
    Query HMDB for metabolite information
    """
    pass


if __name__ == '__main__':
    metabo_data_path = '../data/DyHealthNet/chris_summary_data/metabolites/CHRIS_biocristes7500SumStats.txt'
    metabolite_mapping = read_metabolite_mapping(metabo_data_path)
    hmdb_ids = metabolite_mapping['hmdb_id']
    print(metabolite_mapping.columns)
