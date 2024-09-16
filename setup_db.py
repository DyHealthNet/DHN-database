from settings import *
from cohort_data_format import *
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import DeclarativeBase, sessionmaker
from calculated_edges import add_calculated_edges
from sqlalchemy import URL, text, MetaData, create_engine

from nodes.variants import add_variant_affects_gene, get_cohort_references_variant
from nodes.proteins import read_protein_id_chris, get_protein_nodes, get_protein_interactions
from nodes.variants import get_genomic_variant_nodes, read_rsid_chris, read_variant_meta_file, \
    read_variant_gwas_file
from nodes.metabolites import read_metabolite_mapping, read_hmdb_data, download_metabolite_data, \
    retrieve_assoc_metabolite_nodes
from nodes.phenotypes import download_hpo_ontology, read_hpo_ontology, ontology_data_to_network, snomed_from_hpo, \
    retrieve_disorder_data, retrieve_phenotype_data, get_additional_diseases, get_needed_snomed_ids

from utils.models import *
from utils.query_nedrex import domain_id_to_mondo, get_disorder_data, get_edge_associations, get_gene_data, \
    get_phenotype_data
from utils.database_views import add_views, add_indexes

url = url_object = URL.create(
    "postgresql",
    username=DB_USER,
    password=DB_PASSWORD,
    host=DB_HOST,
    port=DB_PORT,
    database=DB_NAME,
)

engine = create_engine(url)
print("Connected to the database")


def create_tables():
    # does not recreate tables if they already exist
    Base.metadata.create_all(engine)


def delete_tables(session):
    print("Removing all tables from the database.")
    # sql alchemy doesn't support dropping views, so we have to use raw sql
    session.execute(text("DROP MATERIALIZED VIEW IF EXISTS view_description_fts;"))
    session.execute(text("DROP VIEW IF EXISTS view_references_edges;"))
    session.execute(text("DROP VIEW IF EXISTS view_external_nodes;"))
    session.execute(text("DROP MATERIALIZED VIEW IF EXISTS view_associations_edges;"))
    session.commit()
    Base.metadata.drop_all(engine, checkfirst=True)


def add_items(session, items: iter, column: type[DeclarativeBase], filter_args: list, bulk: bool = False):
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
        if exists is not None and DEBUG:
            # remove existing item
            try:
                session.delete(exists)
            except SQLAlchemyError as e:
                session.rollback()
                print("Could not delete existing item: ", e, "moving on...")
        elif exists is not None:
            continue
        try:
            session.add(item)
        except SQLAlchemyError as e:
            session.rollback()
            print("SQLAlchemy Error: ", e)
        except Exception as e:
            session.rollback()
            print("Exception: ", e)
    session.commit()


def add_disorder_data(session, snomed_id_path: str = None, missing_ids: set[str] = None, obs_source: str = None):
    """
    Adds disorder data to the database given a path to a file with snomed ids
    :param session: Database session object
    :param snomed_id_path: str, path to file with snomed ids
    :param obs_source: Describes the source of observations - e.g. CHRIS
    :param missing_ids: Optional - set of omim ids to add to the database. Use this to add missing omim ids from
    i.e. from associations with metabolites
    :return: None
    """
    if missing_ids is None:
        needed_snomed = get_needed_snomed_ids(snomed_id_path)
        data = get_disorder_data(needed_snomed)
        domain_to_mondo = domain_id_to_mondo(data)
    else:
        missing_ids = {f"omim.{x}" for x in missing_ids}
        data = get_disorder_data(missing_ids)
        domain_to_mondo = domain_id_to_mondo(data, 'omim')
        # just to keep downstream code consistent
        needed_snomed = missing_ids

    data = {x['primaryDomainId']: x for x in data}
    xrefs = {mondo: data[mondo]['domainIds'] for mondo in domain_to_mondo.values() if 'domainIds' in data[mondo]}
    display_names = {mondo: data[mondo]['displayName'] for mondo in domain_to_mondo.values()}
    assoc_graph = get_edge_associations(set(domain_to_mondo.values()), edge_type='gene_associated_with_disorder')
    # assoc_graph is filtered for ids that we need, now we can get the data for all genes in the graph since they're
    # all associated with the mondo ids
    gene_info = get_gene_data()
    gene_dict = {x['primaryDomainId']: x for x in gene_info}

    mondo_description = {mondo: data[mondo]['description'] for mondo in domain_to_mondo.values()}

    genes_to_add, disorders, gene_associations, found = retrieve_disorder_data(needed_snomed, domain_to_mondo,
                                                                               mondo_description, xrefs, display_names,
                                                                               gene_dict, assoc_graph, obs_source)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    add_items(session, disorders, Disorder, ['mondo_id'])
    add_items(session, gene_associations, GeneAssocDisorder, ['entrez_id', 'mondo_id'])
    session.commit()
    print(f"Found and successfully added {found} snomed ids with diseases to db")


def add_phenotype_data(session, phenotype_path: str = None, data_dir: str = '../data', missing_ids: list = None,
                       obs_source: str = None):
    """
    Adds phenotype data to the database given a path to a file with phenotype data
    :param obs_source: Describes the source of observations - e.g. CHRIS
    :param session: Database session object
    :param phenotype_path: str, path to file with phenotype data
    :param missing_ids: Optional - set of hpo ids to add to the database. Use this to add missing hpo ids from
    :param data_dir: str, path to the data directory
    :return: None
    """
    # data handling
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    needed_files = [f'{data_dir}/hp.json', f'{data_dir}/phenotype.hpoa']
    if not all([os.path.exists(f) for f in needed_files]):
        download_hpo_ontology(data_dir)

    if not phenotype_path:
        needed_ids = missing_ids
    else:
        needed_ids = get_needed_snomed_ids(phenotype_path)

    hpo_data = read_hpo_ontology(needed_files[0])
    hpo_graph = ontology_data_to_network(hpo_data)

    available_snomed_ids = snomed_from_hpo(hpo_graph, needed_ids)

    available_snomed_ids = {k: v.replace(':', '.').replace('HP', 'hpo') for k, v in available_snomed_ids.items()}
    pheno_data = get_phenotype_data(set(available_snomed_ids.values()))

    additional_data = {item['primaryDomainId']: item for item in pheno_data}
    genes_to_add, phenotypes, disorder_associations, _ = retrieve_phenotype_data(available_snomed_ids, additional_data,
                                                                                 obs_source)

    # since some phenotypes are subtypes of disorders, we only add phenotypes that are
    # not already in the disorder database
    removable_phenotypes = []
    removable_associations = []

    for phenotype in phenotypes:
        snomed = [x for x in phenotype.xrefs if 'snomedct' in x][0]
        items_disorder = session.query(Disorder).filter(Disorder.xrefs.any(snomed)).first()
        if not items_disorder:
            continue
        removable_phenotypes.append(phenotype)

    # also remove associations to phenotypes that are not in the disorder table
    for assoc in disorder_associations:
        if assoc.hpo_id in [x.hpo_id for x in removable_phenotypes]:
            removable_associations.append(assoc)
        if session.query(Disorder).filter_by(mondo_id=assoc.mondo_id).first() is None:
            removable_associations.append(assoc)

    print(f"Removing {len(removable_phenotypes)} phenotypes and {len(removable_associations)} associations that are "
          f"already in the disorder database")
    phenotypes = phenotypes - set(removable_phenotypes)
    disorder_associations = disorder_associations - set(removable_associations)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    add_items(session, phenotypes, Phenotype, ['hpo_id'])
    add_items(session, disorder_associations, DisorderAssocPhenotype, ['mondo_id', 'hpo_id'])

    session.commit()
    print(f"Found and successfully added {len(phenotypes)} snomed ids with phenotypes to db")


def add_cohort_phenotype_data(session, phenotype_path: str = None, obs_source: str = None):
    phenotypes_to_add, phenotype_refs, disorder_refs = cohort_phenotype_data(session, phenotype_path, obs_source)

    add_items(session, phenotypes_to_add, CohortPhenotype, ['cohort_id'])
    add_items(session, phenotype_refs, CohortReferencesPhenotype, ['cohort_id', 'hpo_id'])
    add_items(session, disorder_refs, CohortReferencesDisease, ['cohort_id', 'mondo_id'])
    session.commit()
    print(f"Found and successfully added {len(phenotypes_to_add)} phenotypes from cohort to db")


def add_cohort_metabolite_data(session, metabolite_path: str = None, obs_source: str = None):
    metabolites_to_add, metabolite_refs = cohort_metabolite_data(session, metabolite_path, obs_source)

    add_items(session, metabolites_to_add, CohortMetabolite, ['cohort_id'])
    add_items(session, metabolite_refs, CohortReferencesMetabolite, ['cohort_id', 'hmdb_id'])
    session.commit()
    print(f"Found and successfully added {len(metabolites_to_add)} metabolites from cohort to db")


def add_cohort_protein_data(session, protein_path: str = None, obs_source: str = None):
    proteins_to_add, protein_refs = cohort_protein_data(session, protein_path, obs_source)

    add_items(session, proteins_to_add, CohortProtein, ['cohort_id'])
    add_items(session, protein_refs, CohortReferencesProtein, ['cohort_id', 'uniprot_id'])
    session.commit()
    print(f"Found and successfully added {len(proteins_to_add)} proteins from cohort to db")


def add_cohort_genomic_variants(session, variant_meta_path: str = None, gwas_stats_path: str = None,
                                obs_source: str = None):
    cohort_variants = read_variant_meta_file(variant_meta_path)
    # remove all cohort genomic variants with the same cohort_id
    add_items(session, cohort_variants, CohortVariant, ['cohort_id'], bulk=True)
    session.commit()
    print(f"Added {len(cohort_variants)} cohort genomic variants")
    effect_variant_protein_set, effect_variant_metabolite_set, effect_variant_phenotype_set = read_variant_gwas_file(
        gwas_stats_path)

    def get_cohort_ids(model):
        return {str(row[0]) for row in db_session.query(model.cohort_id).all()}

    rsids_ids = get_cohort_ids(CohortVariant)
    protein_ids = get_cohort_ids(CohortProtein)
    metabolite_ids = get_cohort_ids(CohortMetabolite)
    phenotype_ids = get_cohort_ids(CohortPhenotype)

    variant_protein_filtered = {obj for obj in effect_variant_protein_set if
                                obj.protein_id in protein_ids and obj.variant_id in rsids_ids}
    variant_metabolite_filtered = {obj for obj in effect_variant_metabolite_set if
                                   obj.metabolite_id in metabolite_ids and obj.variant_id in rsids_ids}
    variant_phenotype_filtered = {obj for obj in effect_variant_phenotype_set if
                                  obj.phenotype_id in phenotype_ids and obj.variant_id in rsids_ids}

    add_items(session, variant_metabolite_filtered, EffectVariantMetabolite, ['metabolite_id', 'variant_id'],
              bulk=True)
    add_items(session, variant_phenotype_filtered, EffectVariantPhenotype, ['phenotype_id', 'variant_id'],
              bulk=True)
    add_items(session, variant_protein_filtered, EffectVariantProtein, ['protein_id', 'variant_id'], bulk=True)
    session.commit()

    print(f"Added {len(variant_protein_filtered)} variant-protein associations, "
          f"{len(variant_metabolite_filtered)} variant-metabolite associations, and "
          f"{len(variant_phenotype_filtered)} variant-phenotype associations")

    cohort_references_variant_to_add = get_cohort_references_variant(session, obs_source)
    add_items(session, cohort_references_variant_to_add, CohortReferencesVariant,
              filter_args=["cohort_id", "clinvar_id"])
    session.commit()


def add_protein_data(session, protein_path: str = None, obs_source: str = None, missing_ids: set = None):
    if missing_ids is None:
        protein_ids = read_protein_id_chris(protein_path)
    else:
        protein_ids = missing_ids
    protein_nodes, found_proteins = get_protein_nodes(protein_ids, obs_source)
    needed_ids = {f"uniprot.{uniprot_id}" for uniprot_id in protein_ids}
    print(f"Proteins that couldn't be found: {list(needed_ids - found_proteins)[:5]} and "
          f"{len(needed_ids - found_proteins) - 5} more")

    available_proteins = {x.uniprot_id for x in protein_nodes}
    print(f"Got {len(protein_nodes)} protein nodes")
    protein_interactions = get_protein_interactions(available_proteins)
    add_items(session, protein_nodes, Protein, ['uniprot_id'])
    add_items(session, protein_interactions, ProteinAssocProtein, ['id'])
    session.commit()


def add_metabolite_data(session, metabolite_path: str = None, data_dir: str = '../data', obs_source: str = None):
    hmdb_data_path = f'{data_dir}/hmdb_metabolites.xml'
    download_metabolite_data(data_dir)
    metabolite_mapping = read_metabolite_mapping(metabolite_path)
    unique_metabolites = set()

    # split metabolites that have ; in them
    for metabolite in metabolite_mapping['hmdb_id'].dropna():
        unique_metabolites.update(metabolite.split(';'))

    print(f"Found {len(unique_metabolites)} unique metabolites in the mapping file.")
    omim_diseases = get_additional_diseases(session, obs_source)

    hmdb_mapping = read_hmdb_data(hmdb_data_path, unique_metabolites, omim_ids=omim_diseases,
                                  observation_source=obs_source)
    print(f"Found info for {len(hmdb_mapping)} metabolites in the hmdb data file out of "
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

    print(f"A total of {len(metabolites)} metabolites were found in the mapping file, "
          f"as well as {len(metabolite_protein_associations)} protein associations and "
          f"{len(metabolite_disease_associations)} disease associations.")
    add_items(session, metabolites, Metabolite, ['hmdb_id'])
    add_items(session, metabolite_protein_associations, ProteinAssocMetabolite, ['hmdb_id', 'uniprot_id'])
    add_items(session, metabolite_disease_associations, MetaboliteAssocDisorder, ['hmdb_id', 'mondo_id'])
    session.commit()


def add_genomic_variant_data(session, variant_meta_path: str = None, obs_source: str = None):
    print("Adding genomic variants from NeDRex")
    rs_id_list = read_rsid_chris(variant_meta_path)
    variants_to_add = get_genomic_variant_nodes(rs_id_list, obs_source)
    add_items(session, variants_to_add, Genomic_variant, filter_args=['clinvar_id'])
    session.commit()
    print(f"Added {len(variants_to_add)} genomic variants")

    genomic_variant_ids = {str(row[0]) for row in db_session.query(Genomic_variant.clinvar_id).all()}
    genes_to_add, variant_affects_gene_to_add = add_variant_affects_gene(genomic_variant_ids, obs_source=obs_source)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    session.commit()

    add_items(session, variant_affects_gene_to_add, Variant_affects_gene, filter_args=['entrez_id', 'clinvar_id'])
    session.commit()
    session.commit()
    print("Added variant affects gene edges")


def add_missing(session, data: iter = None, node_type: str = None):
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
        print(f"No missing {node_type} to add to the database.")
        return
    print(f"Got {len(data)} missing {node_type} to add to the database.")
    valid_node_types[node_type](session, missing_ids=data, obs_source='external')


if __name__ == '__main__':
    # Define a session
    Session = sessionmaker(bind=engine)
    db_session = Session()
    delete_tables(db_session)
    dotenv.load_dotenv()
    metadata = MetaData()
    create_tables()
    protein_data_path = PROTEIN_PATH
    pheno_data_path = PHENO_PATH
    metabo_data_path = METABOLITE_PATH
    edges_path = EDGES_PATH
    genomic_variant_meta_path = GENOMIC_VARIANT_META_PATH
    gwas_stats_path = GWAS_STATS_PATH
    data_directory = DATA_DIR

    if not all([pheno_data_path, protein_data_path, metabo_data_path,
                genomic_variant_meta_path, gwas_stats_path, edges_path]):
        raise ValueError("Please provide paths to the phenotype, protein, metabolite, variant and edges files.")

    if not all([os.path.exists(x) for x in [pheno_data_path, protein_data_path, metabo_data_path,
                                            edges_path, genomic_variant_meta_path, gwas_stats_path]]):
        raise ValueError("Some of the provided paths do not exist.")

    print("\nInitialising Layer 2 of database\n")

    add_genomic_variant_data(db_session, genomic_variant_meta_path, obs_source=OBSERVATIONS)
    add_disorder_data(db_session, pheno_data_path, obs_source=OBSERVATIONS)
    add_phenotype_data(db_session, pheno_data_path, obs_source=OBSERVATIONS, data_dir=data_directory)
    add_protein_data(db_session, protein_data_path, obs_source=OBSERVATIONS)
    add_metabolite_data(db_session, metabo_data_path, obs_source=OBSERVATIONS, data_dir=data_directory)

    # second pass for phenotypes
    add_phenotype_data(db_session, pheno_data_path, obs_source='external', data_dir=data_directory)

    # add cohort phenotype data as the mapping is incomplete
    print("\nInitialising Layer 1 of database\n")

    add_cohort_phenotype_data(db_session, pheno_data_path, obs_source=OBSERVATIONS)
    add_cohort_metabolite_data(db_session, metabo_data_path, obs_source=OBSERVATIONS)
    add_cohort_protein_data(db_session, protein_data_path, obs_source=OBSERVATIONS)
    add_cohort_genomic_variants(db_session, genomic_variant_meta_path, gwas_stats_path, obs_source=OBSERVATIONS)

    # add the edges calculated from the available data
    add_calculated_edges(db_session, edges_path, pheno_data_path, protein_data_path, metabo_data_path)

    # count the number of entries in the database
    metadata.reflect(bind=engine)

    # add remaining things (indexes, views)
    add_views(db_session)
    add_indexes(db_session, engine, metadata)

    db_session.close()
    print("Database setup complete.")
