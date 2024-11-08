import sys

from utils.settings import *
from nodes.cohort_nodes import *
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import DeclarativeBase, sessionmaker, Session
from edges.calculated_edges import add_calculated_edges
from sqlalchemy import URL, text, MetaData, create_engine

from nodes.proteins import read_protein_id_chris, get_protein_nodes, get_protein_interactions
from nodes.genes import gene_associations
from nodes.variants import get_genomic_variant_nodes, read_rsid_chris, add_variant_affects_gene
from nodes.metabolites import read_metabolite_mapping, read_hmdb_data, download_metabolite_data, \
    retrieve_assoc_metabolite_nodes
from nodes.phenotypes import download_hpo_ontology, read_hpo_ontology, ontology_data_to_network, \
    retrieve_disorder_data, retrieve_phenotype_data, get_additional_diseases, get_needed_pheno_ids, terms_from_hpo, \
    pheno_ids_from_hpo

from utils.models import *
from utils.query_nedrex import domain_id_to_mondo, get_disorder_data, get_edge_associations, get_gene_data, \
    get_phenotype_data
from utils.database_views import add_views, add_indexes
from utils.logger import get_logger
from analysis.graph_vis import main as graph_vis
from analysis.analyse_db import main as db_stats

logger = get_logger('main')

url = url_object = URL.create(
    "postgresql",
    username=DB_USER,
    password=DB_PASSWORD,
    host=DB_HOST,
    port=DB_PORT,
    database=DB_NAME,
)

engine = create_engine(url)
logger.info("Connected to the database")


def create_tables():
    # does not recreate tables if they already exist
    Base.metadata.create_all(engine)


def delete_tables(session: Session):
    logger.warning("Removing all tables from the database.")
    # sql alchemy doesn't support dropping views, so we have to use raw sql
    session.execute(text("DROP MATERIALIZED VIEW IF EXISTS view_description_fts;"))
    session.execute(text("DROP VIEW IF EXISTS view_references_edges;"))
    session.execute(text("DROP VIEW IF EXISTS view_external_nodes;"))
    session.execute(text("DROP MATERIALIZED VIEW IF EXISTS view_associations_edges;"))
    session.commit()
    Base.metadata.drop_all(engine, checkfirst=True)


def add_items(session: Session, items: iter, column: type[DeclarativeBase], filter_args: list, bulk: bool = False):
    """
    Adds items to the database if they do not already exist
    :param session: Session object
    :param items: iterable of items to add
    :param column: the type of item to add, must be a class from models.py
    :param filter_args: the attributes to filter by to check if the item already exists
    :param bulk: if True, will use bulk_save_objects instead of adding items one by one but this won't check for
    duplicates
    :return: None
    """
    if bulk:
        session.bulk_save_objects(items)
        session.commit()
        return
    for item in items:
        filter_values = {key: getattr(item, key) for key in filter_args}
        exists = session.query(column).filter_by(**filter_values).first()
        if exists is not None and False:#DEBUG:
            # remove existing item
            try:
                session.delete(exists)
            except SQLAlchemyError as e:
                session.rollback()
                logger.warning(f"Could not delete existing item: {e} moving on...")
        elif exists is not None:
            continue
        try:
            session.add(item)
        except SQLAlchemyError as e:
            session.rollback()
            logger.error(f"SQLAlchemy Error: {e}")
        except Exception as e:
            session.rollback()
            logger.error(f"Exception: {e}")
    session.commit()


def add_disorder_data(session: Session, file_path: str = None, missing_ids: set[str] = None, obs_source: str = None):
    """
    Adds disorder data to the database given a path to a file with phenotype reference ids
    :param file_path: str, path to file with phenotype reference ids
    :param session: Database session object
    :param obs_source: Describes the source of observations - e.g. CHRIS
    :param missing_ids: Optional - set of omim ids to add to the database. Use this to add missing omim ids from
    i.e. from associations with metabolites
    :return: None
    """
    if missing_ids is None:
        needed_pheno_ids_no_prefix = get_needed_pheno_ids(file_path)
        needed_pheno_ids = {f"{ID_PREFIX}.{x}" for x in needed_pheno_ids_no_prefix}
        data = get_disorder_data(needed_pheno_ids)
        domain_to_mondo = domain_id_to_mondo(data)
    else:
        needed_pheno_ids_no_prefix = missing_ids
        missing_ids = {f"omim.{x}" for x in missing_ids}
        data = get_disorder_data(missing_ids)
        domain_to_mondo = domain_id_to_mondo(data, 'omim')
        # just to keep downstream code consistent -> not a problem anymore if it was about the snomed related name
        needed_pheno_ids = missing_ids

    data = {x['primaryDomainId']: x for x in data}
    xrefs = {mondo: data[mondo]['domainIds'] for mondo in domain_to_mondo.values() if 'domainIds' in data[mondo]}
    display_names = {mondo: data[mondo]['displayName'] for mondo in domain_to_mondo.values()}
    assoc_graph = get_edge_associations(set(domain_to_mondo.values()), edge_type='gene_associated_with_disorder')
    # assoc_graph is filtered for ids that we need, now we can get the data for all genes in the graph since they're
    # all associated with the mondo ids
    gene_info = get_gene_data()
    gene_dict = {x['primaryDomainId']: x for x in gene_info}

    mondo_description = {mondo: data[mondo]['description'] for mondo in domain_to_mondo.values()}

    genes_to_add, disorders, gene_assocs, found = retrieve_disorder_data(needed_pheno_ids, domain_to_mondo,
                                                                         mondo_description, xrefs, display_names,
                                                                         gene_dict, assoc_graph, obs_source)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    add_items(session, disorders, Disorder, ['mondo_id'])
    add_items(session, gene_assocs, GeneAssocDisorder, ['entrez_id', 'mondo_id'])
    session.commit()
    logger.info(f"Found and successfully added {len(found)} {INPUT_ID_DB} ids with diseases to db")
    # Search for the remaining missing IDs (if any) in NeDRex for phenotype nodes
    curr_missing_ids = list(needed_pheno_ids_no_prefix - found)
    if len(curr_missing_ids) > 0 and missing_ids is None:
        logger.info(f"Adding phenotype...")
        add_phenotype_data(db_session, file_path=None, data_dir=DATA_DIR, missing_ids=curr_missing_ids, obs_source=OBSERVATIONS)

def add_phenotype_data(session: Session, file_path: str = None, data_dir: str = '../data', missing_ids: list = None,
                       obs_source: str = None):
    """
    Adds phenotype data to the database given a path to a file with phenotype data
    :param obs_source: Describes the source of observations - e.g. CHRIS
    :param session: Database session object
    :param file_path: str, path to file with phenotype data if None only missing_ids will be added
    :param missing_ids: Optional - set of phenotype ids to add to the database. Use this to add missing hpo ids from
    :param data_dir: str, path to the data directory
    :return: None
    """

    if HPO_ID_PREFIX is None:
        logger.info("The chose INPUT_DB_ID Type is not present in the HPO API. Therefore no Phenotypes and "
                    "only Disorders will be added.")
        return
    # data handling
    try:
        hpo_graph = terms_from_hpo()
    except Exception:
        logger.info("Could not connect to the HPO API, trying to use flat files")
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        needed_files = [f'{data_dir}/hp.json', f'{data_dir}/phenotype.hpoa']
        if not all([os.path.exists(f) for f in needed_files]):
            download_hpo_ontology(data_dir)

        hpo_data = read_hpo_ontology(needed_files[0])
        hpo_graph = ontology_data_to_network(hpo_data)

    if not file_path:
        needed_ids = missing_ids
    else:
        needed_ids = get_needed_pheno_ids(file_path)

    available_pheno_ids = pheno_ids_from_hpo(hpo_graph, needed_ids)

    available_pheno_ids = {key: [v.replace(':', '.').replace('HP', 'hpo') for v in values]
                           for key, values in available_pheno_ids.items()}
    pheno_data = get_phenotype_data(set(el for v in available_pheno_ids.values() for el in v))

    additional_data = {item['primaryDomainId']: item for item in pheno_data}
    genes_to_add, phenotypes, disorder_associations, _ = retrieve_phenotype_data(available_pheno_ids, additional_data,
                                                                                 obs_source)

    # since some phenotypes are subtypes of disorders, we only add phenotypes that are
    # not already in the disorder database
    removable_phenotypes = set()  #TODO why not use set() directly?
    removable_associations = set()

    if missing_ids is None:
        for phenotype in phenotypes:
            pheno_ids = [x for x in phenotype.xrefs if ID_PREFIX in x][0]
            items_disorder = session.query(Disorder).filter(Disorder.xrefs.any(pheno_ids)).first()
            if not items_disorder:
                continue
            removable_phenotypes.add(phenotype)

    # also remove associations to phenotypes that are not in the disorder table #TODO why not add missing Disorders?
    for assoc in disorder_associations:
        if session.query(Disorder).filter_by(mondo_id=assoc.mondo_id).first() is None:
            removable_associations.add(assoc)
        elif assoc.hpo_id in [x.hpo_id for x in removable_phenotypes]: #TODO faster with static list?
            removable_associations.add(assoc)

    logger.debug(f"Removing {len(removable_phenotypes)} phenotypes and {len(removable_associations)} associations that "
                 f"are already in the disorder database")
    phenotypes = phenotypes - removable_phenotypes
    disorder_associations = disorder_associations - removable_associations

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    add_items(session, phenotypes, Phenotype, ['hpo_id'])
    add_items(session, disorder_associations, DisorderAssocPhenotype, ['mondo_id', 'hpo_id'])

    session.commit()
    logger.info(f"Found and successfully added {len(phenotypes)} {INPUT_ID_DB} ids with phenotypes to db")


def add_cohort_phenotype_data(session: Session, data_path: str = None, obs_source: str = None):
    phenotypes_to_add, phenotype_refs, disorder_refs = cohort_phenotype_data(session, data_path, obs_source)

    add_items(session, phenotypes_to_add, CohortPhenotype, ['cohort_id'])
    add_items(session, phenotype_refs, CohortReferencesPhenotype, ['cohort_id', 'hpo_id'])
    add_items(session, disorder_refs, CohortReferencesDisease, ['cohort_id', 'mondo_id'])
    session.commit()
    logger.info(f"Found and successfully added {len(phenotypes_to_add)} phenotypes from cohort to db")


def add_cohort_metabolite_data(session: Session, data_path: str = None, obs_source: str = None):
    metabolites_to_add, metabolite_refs = cohort_metabolite_data(session, data_path, obs_source)

    add_items(session, metabolites_to_add, CohortMetabolite, ['cohort_id'])
    add_items(session, metabolite_refs, CohortReferencesMetabolite, ['cohort_id', 'hmdb_id'])
    session.commit()
    logger.info(f"Found and successfully added {len(metabolites_to_add)} metabolites from cohort to db")


def add_cohort_protein_data(session: Session, data_path: str = None, obs_source: str = None):
    proteins_to_add, protein_refs = cohort_protein_data(session, data_path, obs_source)

    add_items(session, proteins_to_add, CohortProtein, ['cohort_id'])
    add_items(session, protein_refs, CohortReferencesProtein, ['cohort_id', 'uniprot_id'])
    session.commit()
    logger.info(f"Found and successfully added {len(proteins_to_add)} proteins from cohort to db")


def add_cohort_variants(session: Session, variant_meta_path: str = None, obs_source: str = None):
    cohort_variants, variant_refs = cohort_variant_data(session, variant_meta_path, obs_source)

    add_items(session, cohort_variants, CohortVariant, ['cohort_id'], bulk=True)
    add_items(session, variant_refs, CohortReferencesVariant, filter_args=["cohort_id", "clinvar_id"])
    session.commit()
    logger.info(f"Found and successfully added {len(cohort_variants)} genomic variants from cohort to db")


def add_protein_data(session: Session, file_path: str = None, obs_source: str = None, missing_ids: set = None):
    if missing_ids is None:
        protein_ids = read_protein_id_chris(file_path)
    else:
        protein_ids = missing_ids
    protein_nodes, protein_gene_map = get_protein_nodes(protein_ids, obs_source)

    # get genes not yet in the database
    existing = set(session.query(Gene).filter(Gene.display_name.in_({x[1] for x in protein_gene_map})).all())
    missing = {x[1] for x in protein_gene_map} - {x.entrez_id for x in existing}
    additional_genes, associations = gene_associations(protein_gene_map,
                                                       {x.display_name: x.entrez_id for x in existing}, missing)
    logger.debug(f"Adding {len(additional_genes)} missing genes to the database for {len(associations)} associations")

    needed_ids = {f"uniprot.{uniprot_id}" for uniprot_id in protein_ids}
    found_proteins = {x.uniprot_id for x in protein_nodes}
    logger.info(f"Proteins that couldn't be found: {list(needed_ids - found_proteins)[:5]} and "
                f"{len(needed_ids - found_proteins) - 5} more")

    available_proteins = {x.uniprot_id for x in protein_nodes}
    logger.debug(f"Got {len(protein_nodes)} protein nodes")
    protein_interactions = get_protein_interactions(available_proteins)
    add_items(session, additional_genes, Gene, ['entrez_id'])
    add_items(session, protein_nodes, Protein, ['uniprot_id'])
    add_items(session, associations, ProteinAssocGene, ['uniprot_id', 'entrez_id'])
    add_items(session, protein_interactions, ProteinAssocProtein, ['id'])
    session.commit()


def add_metabolite_data(session: Session, file_path: str = None, data_dir: str = '../data', obs_source: str = None):
    hmdb_data_path = f'{data_dir}/hmdb_metabolites.xml'
    download_metabolite_data(data_dir)
    metabolite_mapping = read_metabolite_mapping(file_path)
    unique_metabolites = set()

    # split metabolites that have ; in them
    for metabolite in metabolite_mapping['hmdb_id'].dropna():
        unique_metabolites.update(metabolite.split(';'))

    logger.debug(f"Found {len(unique_metabolites)} unique metabolites in the mapping file.")
    omim_diseases = get_additional_diseases(session, obs_source)

    hmdb_mapping = read_hmdb_data(hmdb_data_path, unique_metabolites, omim_ids=omim_diseases,
                                  observation_source=obs_source)
    logger.info(f"Found info for {len(hmdb_mapping)} metabolites in the hmdb data file out of "
                f"{len(unique_metabolites)} metabolites in the mapping file.")

    omim_diseases = set()
    for metabolite in hmdb_mapping:
        omim_diseases.update(hmdb_mapping[metabolite]['diseases'])
    omim_diseases = {f"omim.{x}" for x in omim_diseases}
    disease_data = get_disorder_data(omim_diseases)
    omim_mapping = domain_id_to_mondo(disease_data, 'omim')

    metabolites = []
    metabolite_protein_associations = []
    metabolite_disease_associations = []

    prots, diseases = retrieve_assoc_metabolite_nodes(hmdb_mapping)
    missing_proteins = {x for x in prots if session.query(Protein).filter_by(uniprot_id=x).first() is None}
    missing_diseases = {x for x in diseases if session.query(Disorder)
    .filter_by(mondo_id=omim_mapping.get(f"omim.{x}", None)).first() is None}

    add_missing(session, missing_diseases, 'disorders')
    add_missing(session, missing_proteins, 'proteins')

    for metabolite in hmdb_mapping:
        metabolite_name = f"hmdb.{metabolite}"

        metabolites.append(Metabolite(hmdb_id=metabolite_name, display_name=hmdb_mapping[metabolite]['display_name'],
                                      description=hmdb_mapping[metabolite]['description'],
                                      synonyms=hmdb_mapping[metabolite]['synonyms'],
                                      xrefs=hmdb_mapping[metabolite]['xrefs'],
                                      observation_source=hmdb_mapping[metabolite]['observation_source']))

        for disease in hmdb_mapping[metabolite]['diseases']:
            if session.query(Disorder).filter_by(mondo_id=omim_mapping.get(f"omim.{disease}", None)).first() is None:
                continue
            metabolite_disease_associations.append(MetaboliteAssocDisorder(hmdb_id=metabolite_name,
                                                                           mondo_id=omim_mapping[f"omim.{disease}"]))

        for protein in hmdb_mapping[metabolite]['proteins']:
            protein = f"uniprot.{protein}"
            if session.query(Protein).filter_by(uniprot_id=protein).first() is None:
                continue
            metabolite_protein_associations.append(ProteinAssocMetabolite(hmdb_id=metabolite_name, uniprot_id=protein))

    logger.info(f"A total of {len(metabolites)} metabolites were found in the mapping file, "
                f"as well as {len(metabolite_protein_associations)} protein associations and "
                f"{len(metabolite_disease_associations)} disease associations.")
    add_items(session, metabolites, Metabolite, ['hmdb_id'])
    add_items(session, metabolite_protein_associations, ProteinAssocMetabolite, ['hmdb_id', 'uniprot_id'])
    add_items(session, metabolite_disease_associations, MetaboliteAssocDisorder, ['hmdb_id', 'mondo_id'])
    session.commit()


def add_genomic_variant_data(session: Session, file_path: str = None, obs_source: str = None):
    rs_id_list = read_rsid_chris(file_path)
    variants_to_add = get_genomic_variant_nodes(rs_id_list, obs_source)
    add_items(session, variants_to_add, GenomicVariant, filter_args=['clinvar_id'])
    session.commit()
    logger.info(f"Added {len(variants_to_add)} genomic variants")

    genomic_variant_ids = {str(row[0]) for row in db_session.query(GenomicVariant.clinvar_id).all()}
    genes_to_add, variant_affects_gene_to_add = add_variant_affects_gene(genomic_variant_ids, obs_source=obs_source)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    session.commit()

    add_items(session, variant_affects_gene_to_add, VariantAssocGene, filter_args=['entrez_id', 'clinvar_id'])
    session.commit()
    session.commit()
    logger.info("Added variant affects gene edges")


def add_missing(session: Session, data: iter = None, node_type: str = None):
    """
    Adds missing data to the database
    :param session: Session object
    :param data: Data to add
    :param node_type: Type of node to add
    :return: None
    """
    valid_node_types = {
        'proteins': add_protein_data,
        'disorders': add_disorder_data,
        'metabolites': add_metabolite_data,
        'phenotypes': add_phenotype_data,
    }
    if node_type not in valid_node_types:
        raise ValueError(f"Invalid node type: {node_type}")
    if len(data) == 0:
        logger.debug(f"No missing {node_type} to add to the database.")
        return
    logger.info(f"Got {len(data)} missing {node_type} to add to the database.")
    valid_node_types[node_type](session, missing_ids=data, obs_source='external')


def add_node_type(data_path: str = None, data_path_2: str = None) -> bool:
    if not data_path:
        return False

    if not os.path.exists(data_path):
        logger.error(f"Path {data_path} does not exist. Not adding.")
        return False

    if data_path_2 and not os.path.exists(data_path_2):
        logger.error(f"Path {data_path_2} does not exist. Not adding.")
        return False

    return True


def add_layer_node(add_node: bool = False, node_name: str = None, function: callable = None, **kwargs):
    if not add_node:
        logger.debug(f"Node type {node_name} not provided. Therefore, not adding.")
        return
    try:
        logger.info(f"Adding {node_name}...")
        function(**kwargs)
    except Exception as e:
        logger.error(f"Error adding {node_name}: {e}")


if __name__ == '__main__':
    # Define a session
    Session = sessionmaker(bind=engine)
    db_session = Session()

    # delete all tables and recreate them
    delete_tables(db_session)
    create_tables()

    if not all([EDGES_PATH, DATA_DIR]):
        logger.error("Please provide paths to the edges file and data directory.")
        sys.exit(1)

    if not all([os.path.exists(x) for x in [EDGES_PATH, DATA_DIR]]):
        logger.error("Please provide valid paths to the data files.")
        sys.exit(1)

    if ID_PREFIX is None:
        logger.error("Please provide a valid database name that matches the IDs in your phenotype meta file")
        sys.exit(1)

    # returns bool if the node type should be added or not
    add_proteins = add_node_type(PROTEIN_PATH)
    add_phenotypes = add_node_type(PHENO_PATH)
    add_metabolites = add_node_type(METABOLITE_PATH)
    add_variants = add_node_type(VARIANT_META_PATH, EXTRA_EDGES)

    logger.info(f"Will add {'proteins, ' if add_proteins else ''}{'phenotypes, ' if add_phenotypes else ''}"
                f"{'metabolites, ' if add_metabolites else ''}{'variants ' if add_variants else ''}to the db.")

    logger.info("Initialising Layer 2 of database\n")

    # The order at which these functions are called is important
    add_layer_node(add_variants, "genomic variants", add_genomic_variant_data, session=db_session,
                   file_path=VARIANT_META_PATH, obs_source=OBSERVATIONS)

    add_layer_node(add_phenotypes, "disorders", add_disorder_data, session=db_session,
                   file_path=PHENO_PATH, obs_source=OBSERVATIONS)

    # add_layer_node(add_phenotypes, "phenotypes", add_phenotype_data, session=db_session,
    #                file_path=PHENO_PATH, obs_source=OBSERVATIONS, data_dir=DATA_DIR)

    add_layer_node(add_proteins, "proteins", add_protein_data, session=db_session,
                   file_path=PROTEIN_PATH, obs_source=OBSERVATIONS)

    add_layer_node(add_metabolites, "metabolites", add_metabolite_data, session=db_session,
                   file_path=METABOLITE_PATH, obs_source=OBSERVATIONS, data_dir=DATA_DIR)

    # second pass for phenotypes
    if add_phenotypes:
        logger.debug("Doing a second pass for phenotypes to add missing phenotypes") #TODO skip phenotype adding here while still adding their disorder associations
        add_phenotype_data(db_session, PHENO_PATH, obs_source='external', data_dir=DATA_DIR)

    logger.info("Initialising Layer 1 of database\n")

    add_layer_node(add_phenotypes, "cohort phenotypes", add_cohort_phenotype_data, session=db_session,
                   data_path=PHENO_PATH, obs_source=OBSERVATIONS)

    add_layer_node(add_metabolites, "cohort metabolites", add_cohort_metabolite_data, session=db_session,
                   data_path=METABOLITE_PATH, obs_source=OBSERVATIONS)

    add_layer_node(add_proteins, "cohort proteins", add_cohort_protein_data, session=db_session,
                   data_path=PROTEIN_PATH, obs_source=OBSERVATIONS)

    add_layer_node(add_variants, "cohort genomic variants", add_cohort_variants, session=db_session,
                   variant_meta_path=VARIANT_META_PATH, obs_source=OBSERVATIONS)

    # add the edges calculated from the available data
    logger.info("Adding calculated edges...")
    add_calculated_edges(db_session, EDGES_PATH, PHENO_PATH, PROTEIN_PATH, METABOLITE_PATH,
                         VARIANT_META_PATH, EXTRA_EDGES)

    # count the number of entries in the database
    metadata = MetaData()
    metadata.reflect(bind=engine)

    # add remaining things (indexes, views)
    logger.info("Adding views and indexes...")
    add_views(db_session)
    add_indexes(db_session, engine, metadata)

    if VISUALIZE:
        logger.debug("Generating graph visualization...")
        os.makedirs(f'{DATA_DIR}/output') if not os.path.exists(f'{DATA_DIR}/output') else None

        graph_vis(engine, vis_type='html', filename=f'{DATA_DIR}/output/database_graph.html')
        db_stats(db_session, csv_filename=f'{DATA_DIR}/output/database_stats.csv',
                 coverage_filename=f'{DATA_DIR}/output/database_coverage.csv')

    db_session.close()
    logger.info("Database setup complete.")
