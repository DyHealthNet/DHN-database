# load the calculated edges from the file and add them to the database
# steps:
# 1. load the calculated edges from the file
# 2. load info files for label names
# 3. convert label names to proper names
# 4. add the edges to the database

import pandas as pd

from utils.models import *
from sqlalchemy import URL, create_engine, text
from sqlalchemy.orm import sessionmaker
from utils.settings import *

logger = get_logger(__name__)

DB_EDGES = {
    ('protein', 'protein'): (EffectsProteinProtein, 'protein_id_1', 'protein_id_2'),
    ('protein', 'phenotype'): (EffectsProteinPhenotype, 'protein_id', 'phenotype_id'),
    ('phenotype', 'protein'): (EffectsProteinPhenotype, 'phenotype_id', 'protein_id'),
    ('protein', 'metabolite'): (EffectsProteinMetabolite, 'protein_id', 'metabolite_id'),
    ('metabolite', 'protein'): (EffectsProteinMetabolite, 'metabolite_id', 'protein_id'),
    ('protein', 'variant'): (EffectsVariantProtein, 'protein_id', 'variant_id'),

    ('metabolite', 'metabolite'): (EffectsMetaboliteMetabolite, 'metabolite_id_1', 'metabolite_id_2'),
    ('metabolite', 'phenotype'): (EffectsMetabolitePhenotype, 'metabolite_id', 'phenotype_id'),
    ('phenotype', 'metabolite'): (EffectsMetabolitePhenotype, 'phenotype_id', 'metabolite_id'),
    ('metabolite', 'variant'): (EffectsVariantMetabolite, 'metabolite_id', 'variant_id'),

    ('phenotype', 'phenotype'): (EffectsPhenotypePhenotype, 'phenotype_id_1', 'phenotype_id_2'),
    ('phenotype', 'variant'): (EffectsVariantPhenotype, 'phenotype_id', 'variant_id'),

    ('variant', 'metabolite'): (EffectsVariantMetabolite, 'variant_id', 'metabolite_id'),
    ('variant', 'protein'): (EffectsVariantProtein, 'variant_id', 'protein_id'),
    ('variant', 'phenotype'): (EffectsVariantPhenotype, 'variant_id', 'phenotype_id'),
    # ('variant', 'variant'): (EffectsVariantVariant, 'variant_id_1', 'variant_id_2'),
}


def load_files(file_path: str, sep="\t") -> pd.DataFrame | None:
    """
    Load the file from the given path as a pandas DataFrame
    :param file_path: str - path to the file
    :param sep: separator for the file
    :return: Pandas DataFrame
    """
    if not file_path:
        return None
    return pd.read_csv(file_path, sep=sep)


def load_extra(extra_paths: str, sep="\t") -> pd.DataFrame:
    """
    Load the extra files from the given paths as a dictionary of pandas DataFrames
    :param extra_paths: str - path(s), separated by commas
    :param sep: separator for the file
    :return: dict of Pandas DataFrames
    """
    paths = extra_paths.split(',')
    extra_files = []
    for path in paths:
        if not os.path.exists(path):
            logger.error(f"Extra file {path} does not exist, skipping...")
            continue
        extra_files.append(load_files(path, sep))

    assert all([extra_files[0].columns == file.columns for file in extra_files[1:]]), \
        "Extra files do not have the same columns, invalid!"

    return pd.concat(extra_files, ignore_index=True)


def get_labels(label_file: pd.DataFrame, label_type: str = None) -> set[str]:
    """
    Retrieve the set of unique labels from the given file, based on the data type (e.g. protein, metabolite, phenotype)
    The default column names are: protein: protein_id, metabolite: analyte_name, phenotype: label
    :param label_file: Pandas DataFrame
    :param label_type: str - the type of label to retrieve
    :return: set of unique labels
    """
    if not isinstance(label_file, pd.DataFrame):
        return set()
    id_col = COHORT_COLUMNS[label_type]['unique_id']
    # check that the label and id cols exist
    assert id_col in label_file.columns, f"ID column {id_col} not found in the file"
    return set(label_file[id_col].to_list())


def get_xref_rows(session, table: str = 'disorders', xref_col: str = 'omim'):
    """
    Retrieve the rows from the given table that contain the given cross-reference
    :param session: SQLAlchemy session
    :param table: str - database table name
    :param xref_col: str - cross-reference name e.g. omim, snomedct
    :return: list of rows
    """
    # I know this defeats the purpose of SQLAlchemy, but I could not find a way to do this with the ORM
    sql_string = f"""SELECT *
                    FROM {table}
                    WHERE EXISTS (
                        SELECT 1
                        FROM unnest(xrefs) AS xref
                        WHERE xref LIKE '{xref_col}.%'
                    );"""
    return session.execute(text(sql_string)).fetchall()


def filter_exising_ids(session, column: str) -> list:
    """
    Find the existing IDs in the database for the given column
    :param session: SQLAlchemy session
    :param column: str - column name
    :return: list of existing IDs
    """
    return session.query(column).all()


def map_edge(edge: pd.Series, protein_set: set, pheno_set: set, metabo_set: set, variant_set: set) \
        -> tuple[tuple[str, str], tuple[str, str]] | tuple[None, None]:
    """
    Map the source and target of an edge to the appropriate data type given an id and the cohort sets
    :param edge: Pandas Series containing the source and target of the edge
    :param protein_set: set of unique protein IDs from the cohort data
    :param pheno_set: set of unique phenotype labels from the cohort data
    :param metabo_set: set of unique metabolite names from the cohort data
    :return: Tuple containing the mapped source and target, and their respective data types
    """
    maps_and_types = [(protein_set, 'protein'), (pheno_set, 'phenotype'),
                      (metabo_set, 'metabolite'), (variant_set, 'variant')]

    mapped_source = mapped_target = source_type = target_type = None

    for cohort_set, data_type in maps_and_types:
        if not mapped_source:
            if edge['label1'] in cohort_set:
                source_type = data_type
                mapped_source = edge['label1']

        if not mapped_target:
            if edge['label2'] in cohort_set:
                target_type = data_type
                mapped_target = edge['label2']

        if source_type and target_type:
            break

    if not source_type or not target_type:
        return None, None
    return (mapped_source, mapped_target), (source_type, target_type)


def process_chunk(edges_chunk: pd.DataFrame, protein_set: set, phenotype_set: set, metabolite_set: set,
                  variant_set: set) -> tuple[list, list]:
    """
    Process a chunk of edges by mapping and filtering the source and target of the edge, and creating SQLAlchemy objects
    that represent the edge to add to the database.
    Also returns the list of edge types that were added to the database
    :param edges_chunk: Pandas DataFrame containing the chunk of edges
    :param protein_set: set of unique protein IDs from the cohort data
    :param phenotype_set: set of unique phenotype labels from the cohort data
    :param metabolite_set: set of unique metabolite names from the cohort data
    :param variant_set: set of unique variant IDs from the cohort data
    :return: Tuple containing the list of formatted edges and the list of edge types
    """

    def map_and_filter(edge):
        mapped, types = map_edge(edge, protein_set, phenotype_set, metabolite_set, variant_set)
        return mapped, types if types else None

    # Apply the map_and_filter function to the chunk
    mapped_results = edges_chunk.apply(map_and_filter, axis=1)

    # Filter the valid results
    valid_results = [result for result in mapped_results if result[1] is not None]

    if not valid_results:
        return [], []

    # Separate the mapped edges and types
    mapped_edges, types_list = zip(*valid_results)

    # Extract valid rows based on indices of valid results
    valid_indices = [i for i, result in enumerate(mapped_results) if result[1] is not None]
    valid_edges = edges_chunk.iloc[valid_indices]

    # Create a new DataFrame with the mapped data
    formatted_edges = pd.DataFrame({
        'source': [source for source, _ in mapped_edges],
        'target': [target for _, target in mapped_edges],
        'source_type': [source_type for source_type, _ in types_list],
        'target_type': [target_type for _, target_type in types_list],
        'p_value': valid_edges['pval'].values,
        'adjusted_p_value': valid_edges['adj_pval'].values,
        'effect_size': valid_edges['effsize'].values,
        'effect_size_type': valid_edges['effsize_type'].values,
        'test_statistic': valid_edges['test'].values
    })

    formatted_edges_list = []
    for idx, row in formatted_edges.iterrows():
        edge_type = (row['source_type'], row['target_type'])
        if edge_type in DB_EDGES:
            edge_class, source_col, target_col = DB_EDGES[edge_type]
            edge = edge_class(**{
                source_col: row['source'],
                target_col: row['target'],
                'p_value': row['p_value'],
                'adjusted_p_value': row['adjusted_p_value'],
                'effect_size': row['effect_size'],
                'effect_size_type': row['effect_size_type'],
                'test_statistic': row['test_statistic']
            })
            formatted_edges_list.append(edge)

    return formatted_edges_list, [edge.__class__ for edge in formatted_edges_list]


def format_edges(session, edges: pd.DataFrame, protein_set: set, phenotype_set: set, metabolite_set: set,
                 variant_set: set) -> None:
    """
    Format the edges and add them to the database in chunks. Deletes the formatted edges after adding them to the
    database to save memory. The chunk size can be adjusted in the settings.
    :param session: SQLAlchemy session
    :param edges: Pandas DataFrame containing the edges
    :param protein_set: set of unique protein IDs from the cohort data
    :param phenotype_set: set of unique phenotype labels from the cohort data
    :param metabolite_set: set of unique metabolite names from the cohort data
    :param variant_set: set of unique variant IDs from the cohort data
    :return: None
    """
    all_edge_types = {edge_type[0]: 0 for edge_type in DB_EDGES.values()}
    num_edge_types = {}
    chunk_size = CHUNK_SIZE
    num_chunks = (len(edges) // chunk_size) + 1

    if DEBUG:
        chunk_size = 10_000
        num_chunks = 1
        edges = edges.sample(frac=1)

    for i in range(num_chunks):
        start_index = i * chunk_size
        end_index = min((i + 1) * chunk_size, len(edges))
        edges_chunk = edges.iloc[start_index:end_index]

        # Process the current chunk
        formatted_edges_list, edge_types = process_chunk(edges_chunk, protein_set, phenotype_set,
                                                         metabolite_set, variant_set)

        add_success = add_edges(session, formatted_edges_list)
        del formatted_edges_list
        if not add_success:
            logger.error("There was a problem adding the edges to the database, exiting...")
            return
        logger.info(f"Chunk {i + 1}/{num_chunks} added successfully")
        # do the value counts of the edges and add them to a running total
        chunk_edge_types = pd.Series(edge_types).value_counts().to_dict()
        num_edge_types = {edge_type: num_edge_types_value + chunk_edge_types.get(edge_type, 0)
                          for edge_type, num_edge_types_value in all_edge_types.items()}

    for edge_type, count in num_edge_types.items():
        logger.debug(f"Added {count} edges of type {edge_type.__name__}")
    return


def add_edges(session, edges: list[Base]) -> bool:
    """
    Add the given list of edges to the database in bulk
    :param session: SQLAlchemy session
    :param edges: list of SQLAlchemy objects representing the edges
    :return: bool - True if the edges were added successfully, False otherwise
    """
    try:
        session.bulk_save_objects(edges)
        session.commit()
    except Exception as e:
        session.rollback()
        logger.error(f"A problem occurred while adding edges: {e}")
        return False
    return True


def add_calculated_edges(session, edges_path: str,
                         pheno_data_path: str | None,
                         protein_data_path: str | None,
                         metabo_data_path: str | None,
                         variant_data_path: str | None,
                         extra_edge_paths: str | None = None) -> None:
    """
    Main function to add the calculated edges to the database.
    :param session: SQLAlchemy session
    :param edges_path: str - path to the calculated edges file
    :param pheno_data_path: str - path to the phenotype data file
    :param protein_data_path: str - path to the protein data file
    :param metabo_data_path: str - path to the metabolite data file
    :param variant_data_path: str - path to the variant data file
    :param extra_edge_paths: str - path to the extra data files (e.g. variant edges)
    :return: None
    """
    phenotypes = load_files(pheno_data_path)
    proteins = load_files(protein_data_path)
    metabolites = load_files(metabo_data_path)
    variants = load_files(variant_data_path)

    edges = load_files(edges_path, sep=",")
    edges = edges.iloc[:, 1:]  # removing unnamed column
    if extra_edge_paths:
        logger.debug("Loading extra data files")
        extra_data = load_extra(extra_edge_paths)
        edges = pd.concat([edges, extra_data], ignore_index=True)

    # get base labels
    pheno_set = get_labels(phenotypes, 'phenotype')
    protein_set = get_labels(proteins, 'protein')
    metabo_set = get_labels(metabolites, 'metabolite')
    variant_set = get_labels(variants, 'variant')
    logger.debug("All cohort sets loaded successfully")

    format_edges(session, edges, protein_set, pheno_set, metabo_set, variant_set)
    # formatted_edges = format_edges(edges[edges['pval'] <= 0.05], protein_set, pheno_map, metabo_map, disorder_map)


if __name__ == '__main__':
    url = url_object = URL.create(
        "postgresql",
        username=DB_USER,
        password=DB_PASSWORD,
        host=DB_HOST,
        port=DB_PORT,
        database=DB_NAME
    )
    engine = create_engine(url)
    Session = sessionmaker(bind=engine)
    db_session = Session()

    edges_path = '../data/scores_including_tests.csv'
    pheno_data_path = '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv'
    protein_data_path = None
    metabo_data_path = '../data/DyHealthNet/chris_summary_data/metabolites/CHRIS_biocristes7500SumStats.txt'
    variant_data_path = None

    add_calculated_edges(db_session, edges_path, pheno_data_path, protein_data_path,
                         metabo_data_path, variant_data_path)
