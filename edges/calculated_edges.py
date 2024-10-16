# load the calculated edges from the file and add them to the database
# steps:
# 1. load the calculated edges from the file
# 2. load info files for label names
# 3. convert label names to proper names
# 4. add the edges to the database
from io import StringIO

import pandas as pd

from utils.models import *
from sqlalchemy import URL, create_engine, text
from sqlalchemy.orm import sessionmaker, Session
from utils.settings import *

logger = get_logger(__name__)


DB_EDGES = {
    ('phenotype', 'variant'): "effects_variant_phenotype",
    ('variant', 'phenotype'): "effects_variant_phenotype",
    ('variant', 'metabolite'): "effects_variant_metabolite",
    ('metabolite', 'variant'): "effects_variant_metabolite",
    ('variant', 'protein'): "effects_variant_protein",
    ('protein', 'variant'): "effects_variant_protein",

    ('protein', 'protein'): "effects_protein_protein",
    ('protein', 'phenotype'): "effects_protein_phenotype",
    ('phenotype', 'protein'): "effects_protein_phenotype",
    ('protein', 'metabolite'): "effects_protein_metabolite",
    ('metabolite', 'protein'): "effects_protein_metabolite",

    ('metabolite', 'metabolite'): "effects_metabolite_metabolite",
    ('metabolite', 'phenotype'): "effects_metabolite_phenotype",
    ('phenotype', 'metabolite'): "effects_metabolite_phenotype",

    ('phenotype', 'phenotype'): "effects_phenotype_phenotype",
}

# this gives us information which data type comes first in the column order
EDGE_ORDER = {'variant': 3, 'protein': 2, 'metabolite': 1, 'phenotype': 0}


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


def load_extra(extra_paths: str, sep="\t") -> StringIO:
    """
    Load the extra files from the given paths as a dictionary of pandas DataFrames
    :param extra_paths: str - path(s), separated by commas
    :param sep: separator for the file
    :return: dict of Pandas DataFrames
    """
    paths = extra_paths.split(',')
    buffer = StringIO()

    for path in paths:
        if not os.path.exists(path):
            logger.error(f"Extra file {path} does not exist, skipping...")
            continue
        raw_df = load_files(path, sep)
        try:
            raw_df = raw_df[['label1', 'label2', 'pval', 'adj_pval', 'effsize', 'effsize_type', 'test']]
        except KeyError:
            logger.error(f"Extra file {path} is missing columns, skipping...")
            continue

        # replace column pval and effsize with an escaped version
        raw_df['pval'] = raw_df['pval'].apply(lambda x: r'"{\"full\": %s}"' % x)
        raw_df['effsize'] = raw_df['effsize'].apply(lambda x: r'"{\"full\": %s}"' % x)
        raw_df.to_csv(buffer, sep=',', index=True, header=False, quoting=3, lineterminator='\n')
        # move the buffer to the end of the file
        buffer.seek(0, 2)

    return buffer


def swap_labels(s, label1, label2):
    index1 = s.find(label1)
    index2 = s.find(label2)

    if index1 == -1 or index2 == -1:
        return s

    before = s[:min(index1, index2)]
    between = s[min(index1 + len(label1), index2 + len(label2)):max(index1, index2)]
    after = s[max(index1 + len(label1), index2 + len(label2)):]

    if index1 < index2:
        return before + label2 + between + label1 + after
    return before + label1 + between + label2 + after


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


def map_edge(edge: tuple, protein_set: set, pheno_set: set, metabo_set: set, variant_set: set) \
        -> tuple[str, str, bool] | None:
    """
    Map the source and target of an edge to the appropriate data type given an id and the cohort sets
    :param variant_set: set of unique variant IDs from the cohort data
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
            if edge[0] in cohort_set:
                source_type = data_type
                mapped_source = edge[0]

        if not mapped_target:
            if edge[1] in cohort_set:
                target_type = data_type
                mapped_target = edge[1]

        if source_type and target_type:
            break

    if not source_type or not target_type:
        return None

    swap = False
    if EDGE_ORDER[source_type] < EDGE_ORDER[target_type]:
        swap = True
    return source_type, target_type, swap


def process_file(edges, protein_set: set, phenotype_set: set, metabolite_set: set,
                 variant_set: set) -> dict:
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
    all_edge_types = {edge_type: StringIO() for edge_type in DB_EDGES.values()}

    def map_and_filter(edge: tuple) -> tuple[tuple[str, str], tuple[str, str], bool] | tuple[None, None, bool]:
        mapped, types, swap = map_edge(edge, protein_set, phenotype_set, metabolite_set, variant_set)
        return mapped, types if types else None, swap

    edges.readline()
    for line in edges.readlines():
        if "nan" in line:
            continue
        line_split = line.split(',')
        source, dest = line_split[1], line_split[2]
        source_map, dest_map, swap = map_and_filter((source, dest))
        if source_map is None or dest_map is None:
            continue
        edge_map = (source_map, dest_map)
        # swap the labels to match the order in the database
        if swap:
            line = swap_labels(line, source, dest)
        all_edge_types[DB_EDGES[edge_map]].write(line)

    logger.debug("Finished processing edges")
    return all_edge_types


def format_edges(session: Session, edges, protein_set: set, phenotype_set: set, metabolite_set: set,
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

    formatted_edges = process_file(edges, protein_set, phenotype_set, metabolite_set, variant_set)

    add_success = add_edges(session, formatted_edges)
    del formatted_edges
    if not add_success:
        logger.error("There was a problem adding the edges to the database, exiting...")
        return
    logger.info(f"File added successfully")

    # for edge_type, count in num_edge_types.items():
    #    logger.debug(f"Added {count} edges of type {edge_type.__name__}")
    return


def copy_from_buffer(session, edge_type, edge_file):
    # Get the raw psycopg2 connection from SQLAlchemy
    raw_connection = session.connection().connection
    with raw_connection.cursor() as cursor:
        edge_file.seek(0)
        copy_sql = f"COPY {edge_type} FROM STDIN WITH CSV HEADER QUOTE '\"' DELIMITER ',' ESCAPE '\\'"
        cursor.copy_expert(copy_sql, edge_file)
        raw_connection.commit()


def count_rows(buffer):
    buffer.seek(0)
    row_count = buffer.getvalue().count('\n')
    return row_count


def add_edges(session: Session, edges: dict) -> bool:
    """
    Add the given list of edges to the database in bulk
    :param session: SQLAlchemy session
    :param edges: list of SQLAlchemy objects representing the edges
    :return: bool - True if the edges were added successfully, False otherwise
    """
    try:
        if DEBUG:
            # clear the tables
            for edge_type in DB_EDGES.values():
                session.execute(text(f"TRUNCATE TABLE {edge_type}"))
        for edge_type, edge_file in edges.items():
            if edge_file is None:
                continue
            copy_from_buffer(session, edge_type, edge_file)
            edge_count = count_rows(edge_file)
            if edge_count > 0:
                logger.debug(f"Finished adding {edge_count} {edge_type} edges")
    except Exception as e:
        session.rollback()
        logger.error(f"A problem occurred while adding edges: {e}")
        return False
    return True


def add_calculated_edges(session: Session,
                         edges_path: str,
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

    # get base labels
    pheno_set = get_labels(phenotypes, 'phenotype')
    protein_set = get_labels(proteins, 'protein')
    metabo_set = get_labels(metabolites, 'metabolite')
    variant_set = get_labels(variants, 'variant')
    logger.debug("All cohort sets loaded successfully")

    edges = open(edges_path, 'r')
    format_edges(session, edges, protein_set, pheno_set, metabo_set, variant_set)

    if extra_edge_paths:
        logger.debug("Loading extra data files")
        extra_data = load_extra(extra_edge_paths)
        extra_data.seek(0)
        format_edges(session, extra_data, protein_set, pheno_set, metabo_set, variant_set)


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

    edges_path = EDGES_PATH
    pheno_data_path = None
    protein_data_path = "../" + PROTEIN_PATH
    metabo_data_path = METABOLITE_PATH
    variant_data_path = None

    add_calculated_edges(db_session, edges_path, pheno_data_path, protein_data_path,
                         metabo_data_path, variant_data_path)
