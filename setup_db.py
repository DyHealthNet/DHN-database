import networkx as nx
from settings import *
from nedrex.core import iter_nodes, iter_edges, get_node_types
from collections import Counter
from testcases import *
from cohort_data_format import *
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import DeclarativeBase
from calculated_edges import add_calculated_edges
from sqlalchemy import URL, text, Table, MetaData, distinct
from protein_mapping import read_proteinID_chris, get_protein_nodes
from genomic_variant_to_clinvar import get_genomic_variant_nodes, read_rsid_chris, read_variant_meta_file, \
    read_variant_gwas_file
from metabolite_mapping import read_metabolite_mapping, read_hmdb_data, download_metabolite_data, \
    retrieve_assoc_metabolite_nodes
from hpo_mapping import download_hpo_ontology, read_hpo_ontology, ontology_data_to_network, snomed_from_hpo
from query_nedrex import get_needed_snomed_ids, domain_id_to_mondo, get_disorder_data, get_edge_associations, \
    get_harmonizome_data, get_gene_data, get_phenotype_data
from database_views import add_views, add_indexes
import csv

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
    session.execute(text("DROP VIEW IF EXISTS external_node_ids;"))
    session.execute(text("DROP MATERIALIZED VIEW IF EXISTS view_associations_edges;"))
    session.commit()
    Base.metadata.drop_all(engine, checkfirst=True)


def mondo_in_association_graph(mondo_id: str, assoc_graph: nx.Graph) -> tuple[list[str], list[str]] | tuple[None, None]:
    """
    Retrieves the genes associated with a mondo id from the association graph
    :param mondo_id: A mondo id, can be None
    :param assoc_graph: The association graph from NEDRex
    :return: genes associated with the mondo id, source of the data, harmonizome data if not in the association graph
    """
    if mondo_id is None:
        return None, None
    if mondo_id in assoc_graph:
        genes = []
        sources = []
        for edge in assoc_graph.edges(mondo_id, data=True):
            genes.append(edge[1])
            sources.append(edge[2]['source'][0])
    else:
        harm = get_harmonizome_data(mondo_id)
        if harm is None:
            return None, None
        genes = []
        sources = []
        for key, value in harm.items():
            genes.extend(value)
            sources.extend([key] * len(value))
    return genes, sources


def retrieve_disorder_data(needed_snomed: set[str], snomed_to_mondo: dict[str, str], descriptions: dict[str, str],
                           xrefs: dict, display_names: dict, gene_info: dict, assoc_graph: nx.Graph,
                           obs_source: str = None) -> tuple[set, set, set, int]:
    """
    Queries the needed snomed ids and retrieves the associated genes and disorders from NeDRex
    :param display_names: Display names for the mondo ids
    :param gene_info: Information about the genes needed for the database (display name, synonyms, etc.)
    :param xrefs: cross references for the mondo ids to other databases
    :param descriptions: descriptions for the mondo ids
    :param needed_snomed: snomed ids in the dataset
    :param snomed_to_mondo: map from snomed to mondo ids from nedrex
    :param assoc_graph: association graph from nedrex of mondo ids to genes
    :param obs_source: Describes the source of observations - e.g. CHRIS
    :return: list of genes to add, list of disorders to add, list of gene associations to add,
    number of snomed ids found
    """
    gene_associations = set()
    genes_to_add = set()
    disorders = set()
    found = 0
    for snomed in needed_snomed:
        # check if ; in snomed id and if so, do this for all ids
        snomed_ids = str(snomed).split(';')
        for snomed_id in snomed_ids:
            mondo_id = snomed_to_mondo.get(snomed_id)
            xref = xrefs.get(mondo_id, None)
            description = descriptions.get(mondo_id, None)
            display_name = display_names.get(mondo_id, None)
            genes, sources = mondo_in_association_graph(mondo_id, assoc_graph)
            if genes is None:
                continue
            # add genes to set
            for gene in genes:
                gene_data = gene_info.get(gene, None)
                if gene_data is None:
                    new_gene = Gene(entrez_id=gene, observation_source='external')
                    genes_to_add.add(new_gene)
                    continue
                new_gene = Gene(entrez_id=gene_data['primaryDomainId'],
                                display_name=gene_data['displayName'],
                                description=gene_data['description'],
                                synonyms=gene_data['synonyms'],
                                chromosome=gene_data['chromosome'],
                                observation_source='external')
                genes_to_add.add(new_gene)

            disorders.add(
                Disorder(mondo_id=mondo_id,
                         display_name=display_name,
                         xrefs=xref,
                         description=description,
                         observation_source=obs_source))
            # add gene associations to set for each source
            gene_associations.update([GeneAssocDisorder(entrez_id=gene, mondo_id=mondo_id, edge_source=source)
                                      for gene, source in zip(genes, sources)])
            found += 1
        continue
    return genes_to_add, disorders, gene_associations, found


def retrieve_phenotype_data(available_ids: dict, additional_data: dict, obs_source: str = None) \
        -> tuple[set, set, set, int]:
    """
    Retrieve the phenotype data for the needed snomed ids
    # HPO conversion: Pathway
    # HPO data (HPO_ID ----> SNOMED_ID) - look for needed SNOMED IDs
    # -> map to Mondo (SNOMED_ID -- OMIM_ID/ORPHA_ID --> Mondo_ID)
    #

    :param obs_source: Describes the source of observations - e.g. CHRIS
    :param additional_data: dictionary with additional data for the hpo ids, must be a dictionary with hpo ids as keys
    :param available_ids: dictionary with snomed ids as keys and hpo ids as values
    :return: dictionary with the phenotype data
    """
    found = 0
    genes_to_add = set()
    phenotypes = set()
    disorder_associations = set()

    available_snomed_ids = available_ids
    # go through all the nodes in the HPO graph and find the ones that have xrefs to SNOMED

    print(f'Found {len(available_snomed_ids)} snomed ids in the HPO ontology')

    # get edge associations for disorder_has_phenotype
    assoc_graph = get_edge_associations(set(available_snomed_ids.values()), edge_type='disorder_has_phenotype')

    # find the genes that are associated with the mondo ids
    for snomed_id, hpo_id in available_snomed_ids.items():
        snomed_id = f"snomedct.{snomed_id}"
        if additional_data.get(hpo_id, None) is None:
            continue
        phenotype_data = additional_data[hpo_id]
        phenotypes.add(Phenotype(hpo_id=hpo_id,
                                 xrefs=set(phenotype_data['domainIds'] + [snomed_id]),
                                 description=phenotype_data['description'],
                                 synonyms=phenotype_data['synonyms'],
                                 display_name=phenotype_data['displayName'],
                                 observation_source=obs_source))
        if hpo_id not in assoc_graph:
            continue
        # get the disorder ids associated with the hpo id
        for edge in assoc_graph.edges(hpo_id, data=True):
            disorder = edge[1]
            source = edge[2]['source'][0]
            new_assoc = DisorderAssocPhenotype(mondo_id=disorder, hpo_id=hpo_id, edge_source=source)
            disorder_associations.add(new_assoc)

        found += 1
    return genes_to_add, phenotypes, disorder_associations, found


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

    #         hpo_id = hpo_id.replace(':', '.').replace('HP', 'hpo')
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


def add_protein_data(session, proteinData_path = None, obs_source = None, missing_ids = None):
    if missing_ids is None:
        proteinIds = read_proteinID_chris(proteinData_path)
    else:
        proteinIds = missing_ids
    proteinNodes, found_proteins = get_protein_nodes(proteinIds, obs_source)
    needed_ids = {f"uniprot.{uniprot_id}" for uniprot_id in proteinIds}
    print(f"Proteins that couldn't be found: {list(needed_ids - found_proteins)[:5]} and "
          f"{len(needed_ids - found_proteins) - 5} more")

    available_proteins = {x.uniprot_id for x in proteinNodes}
    print(f"Got {len(proteinNodes)} protein nodes")
    proteinInteractions = get_protein_interactions(available_proteins)
    add_items(session, proteinNodes, Protein, ['uniprot_id'])
    add_items(session, proteinInteractions, ProteinAssocProtein, ['id'])
    session.commit()


def get_protein_interactions(proteinIds):
    # prefixed_proteinIds = {f"uniprot.{entry}" for entry in proteinIds}
    # retrieve_interacting_proteins_neo4j(proteinIds)
    assoc_graph = get_edge_associations(proteinIds, edge_type='protein_interacts_with_protein',
                                        direction='undirected')
    proteinInteractions = []
    for edge in assoc_graph.edges():
        uniprot_id_memberOne = edge[0]
        uniprot_id_memberTwo = edge[1]
        if (uniprot_id_memberTwo in proteinIds and uniprot_id_memberOne in proteinIds):
            proteinInteractions.append(ProteinAssocProtein(uniprot_id_1=uniprot_id_memberOne,
                                                           uniprot_id_2=uniprot_id_memberTwo))
    return proteinInteractions


def get_additional_diseases(session, obs_source: str = None):
    # I know this defeats the purpose of SQLAlchemy but I could not find a way to do this with the ORM
    sql_string = f"""SELECT xrefs
                    FROM disorder
                    WHERE EXISTS (
                        SELECT 1
                        FROM unnest(xrefs) AS xref
                        WHERE xref LIKE 'omim.%'
                    ) AND observation_source = '{obs_source}';"""
    return {x for x in session.execute(text(sql_string)).fetchall() for x in x[0] if x.startswith('omim.')}


def add_metabolite_data(session, metabolite_path, data_dir: str = '../data', obs_source: str = None):
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


def add_cohort_genomic_variants(session, genomic_variant_meta_path, gwas_stats_path, obs_source):
    gwas_stats_path = gwas_stats_path
    cohort_variants = read_variant_meta_file(genomic_variant_meta_path)
    # remove all cohort genomic variants with the same cohort_id
    id_counts = Counter(x.cohort_id for x in cohort_variants)
    cohort_variants = [x for x in cohort_variants if id_counts[x.cohort_id] == 1]
    add_items(session, cohort_variants, CohortGenomicVariant, ['cohort_id'], bulk=True)
    session.commit()
    print(f"Added {len(cohort_variants)} cohort genomic variants")
    effectVariantProteinSet, effectVariantMetaboliteSet, effectVariantPhenotypeSet = read_variant_gwas_file(
        gwas_stats_path)

    rsids_ids = {str(row[0]) for row in db_session.query(CohortGenomicVariant.cohort_id).all()}
    phenotype_ids = {str(row[0]) for row in db_session.query(CohortPhenotype.cohort_id).all()}
    protein_ids = {str(row[0]) for row in db_session.query(CohortProtein.cohort_id).all()}
    metabolite_ids = {str(row[0]) for row in db_session.query(CohortMetabolite.cohort_id).all()}
    effectVariantProteinSet_filtered = {obj for obj in effectVariantProteinSet if
                                        obj.protein_id in protein_ids and obj.variant_id in rsids_ids}
    effectVariantMetaboliteSet_filtered = {obj for obj in effectVariantMetaboliteSet if
                                           obj.metabolite_id in metabolite_ids and obj.variant_id in rsids_ids}
    effectVariantPhenotypeSet_filtered = {obj for obj in effectVariantPhenotypeSet if
                                          obj.phenotype_id in phenotype_ids and obj.variant_id in rsids_ids}
    add_items(session, effectVariantMetaboliteSet_filtered, EffectVariantMetabolite, ['metabolite_id', 'variant_id'], bulk=True)
    add_items(session, effectVariantPhenotypeSet_filtered, EffectVariantPhenotype, ['phenotype_id', 'variant_id'], bulk=True)
    add_items(session, effectVariantProteinSet_filtered, EffectVariantProtein, ['protein_id', 'variant_id'], bulk=True)
    session.commit()
    print(f"Added {len(effectVariantProteinSet_filtered)} variant-protein associations, "
          f"{len(effectVariantMetaboliteSet_filtered)} variant-metabolite associations, and "
          f"{len(effectVariantPhenotypeSet_filtered)} variant-phenotype associations")
    cohort_references_variant_to_add = get_cohort_references_variant(session, obs_source)
    add_items(session, cohort_references_variant_to_add, CohortReferencesVariant,
              filter_args=["cohort_id", "clinvar_id"])
    session.commit()


def get_cohort_references_variant(session, obs_source):
    genomic_variants = session.query(Genomic_variant).all()
    newCohortReferencesSet = set()
    #existing_rsids_display = {str(row[0]) for row in session.query(CohortGenomicVariant.display_name).all()}
    existing_cohort_id = {x[0] for x in session.query(CohortGenomicVariant.cohort_id).all()}
    #existing_clinvar_ids = {str(row[0]) for row in session.query(Genomic_variant.clinvar_id).all()}
    for variant in genomic_variants:
        variant_domain_ids = variant.xrefs.replace(",", "").replace("[", "").replace("]", "").replace("'","").split()
        ref_cohort_id =str(variant.chromosome + ":" + variant.position + ":" + variant.referenceSequence + ">" + variant.alternativeSequence)
        dbsnp_id = next((id.replace("dbsnp.", "rs") for id in variant_domain_ids if "dbsnp." in id), None)
        clinvar_id = variant.clinvar_id
        #next((id for id in variant_domain_ids if "clinvar." in id), None)
        if (ref_cohort_id in existing_cohort_id):
            newCohortReferencesVariant = CohortReferencesVariant(
                cohort_id=ref_cohort_id,
                clinvar_id=clinvar_id
            )
            newCohortReferencesSet.add(newCohortReferencesVariant)

    return (newCohortReferencesSet)


def add_genomic_variants_neddrex(session, genomic_variant_meta_path, obs_source):
    print("Adding genomic variants from NeDRex")
    rs_id_list = read_rsid_chris(genomic_variant_meta_path)
    variants_to_add = get_genomic_variant_nodes(rs_id_list, obs_source)
    add_items(session, variants_to_add, Genomic_variant, filter_args=['clinvar_id'])
    session.commit()
    genomic_variant_ids = {str(row[0]) for row in db_session.query(Genomic_variant.clinvar_id).all()}
    add_genomic_variant_edge_variant_affects_gene(db_session, genomic_variant_ids, obs_source=obs_source)
    session.commit()


def add_genomic_variant_edge_variant_affects_gene(session, clinvarIds, obs_source):
    variant_affects_gene_graph = get_edge_associations(node_ids=clinvarIds, edge_type='variant_affects_gene',
                                                       direction='directed')
    variant_affects_gene_dict = {}
    id_list_total = []  # to check if is in database
    variant_edge_list = []
    variant_affects_gene_to_add = set()
    for edge in variant_affects_gene_graph.edges(data=True):
        variant_affects_gene_dict[edge[0]] = edge[1]
        variant_edge_list.append(edge)
        sourceDomainId = edge[0]
        targedDomainId = edge[1]
        id_list_total.append(sourceDomainId)
        id_list_total.append(
            targedDomainId)  # because graph is directed either first or second entry contains the entrez id

        if sourceDomainId.startswith("entrez."):
            entrez_id = sourceDomainId
            variant = targedDomainId
        if targedDomainId.startswith("entrez."):
            entrez_id = targedDomainId
            variant = sourceDomainId
        variant_affects_gene_edge = Variant_affects_gene(clinvar_id=variant,
                                                         entrez_id=entrez_id)
        variant_affects_gene_to_add.add(variant_affects_gene_edge)
    entrez_ids = {id for id in id_list_total if id.startswith("entrez.")}
    genomic_variant_node_generator = iter_nodes('genomic_variant')
    genomic_variant_node_generator = iter_nodes('gene')
    genes_to_add = set()
    for node in genomic_variant_node_generator:
        if (node['primaryDomainId'] in id_list_total):
            new_gene = Gene(entrez_id=node['primaryDomainId'],
                            display_name=node['displayName'],
                            description=node['description'],
                            synonyms=node['synonyms'],
                            chromosome=node['chromosome'],
                            observation_source='external')
            genes_to_add.add(new_gene)
    add_items(session, genes_to_add, Gene, ['entrez_id'])
    session.commit()

    add_items(session, variant_affects_gene_to_add, Variant_affects_gene,
              filter_args=['entrez_id', 'clinvar_id'])
    session.commit()

    # stmt = select([genomic_variants.c.id]).where(genomic_variants.c.id == id_to_check)
    # result = db_session.execute(stmt).first()
    return None


def add_genomic_variants_linking_from_gene(session, entrez_ids: set[str] = None, observation_source: str = None):
    variant_affects_gene_graph = get_edge_associations(node_ids=entrez_ids, edge_type='variant_affects_gene',
                                                       direction='directed')  # Graph with 1511628 nodes and 1539719 edges
    print(f"Got {len(variant_affects_gene_graph)} edge associations for variant_affects_gene.")
    variant_affects_gene_dict = {}
    clinvar_ids = set()
    for edge in variant_affects_gene_graph.edges(data=True):
        clinvar_ids.add(edge[0])
        variant_affects_gene_dict[edge[0]] = edge[1]

    pattern = r'^[^.]*\.'
    variants_to_add = []
    genomic_variant_node_generator = iter_nodes('genomic_variant')
    for node in genomic_variant_node_generator:
        if node['primaryDomainId'] in clinvar_ids:
            newVariant = Genomic_variant(clinvar_id=node['primaryDomainId'],  # linvar.17735
                                         alternativeSequence=node['alternativeSequence'],  # 'T',
                                         chromosome=node['chromosome'],  # 'NW_009646201.1',
                                         dataSources=node['dataSources'],  # ['clinvar'],
                                         xrefs=node['domainIds'],  # ['clinvar.17735', 'dbsnp.1556058284']
                                         position=node['position'],  # 83614,
                                         referenceSequence=node['referenceSequence'],  # 'TC',
                                         type=node['type'],  # 'GenomicVariant'
                                         variantType=node['variantType'])  # 'Deletion'})
            variants_to_add.append(newVariant)

    add_items(session, variants_to_add, Genomic_variant, filter_args=['clinvar_id'])
    variant_affects_gene_to_add = []
    genes = []
    available_variants = {variant.clinvar_id for variant in variants_to_add}
    for variant, gene in variant_affects_gene_dict.items():
        if (variant in available_variants and gene in entrez_ids):
            variant_affects_gene_to_add.append(Variant_affects_gene(clinvar_id=variant, entrez_id=gene))
    add_items(session, variant_affects_gene_to_add, Variant_affects_gene, filter_args=['entrez_id', 'clinvar_id'])
    session.commit()


def add_missing(session, data, node_type):
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


def countEntries(session, metadata):
    table_counts = {}
    # Iterate over each table in the metadata
    for table_name in metadata.tables:
        table = Table(table_name, metadata, autoload_with=engine)
        count = session.query(table).count()
        table_counts[table_name] = count
    for table_name, count in table_counts.items():
        print(f"Table {table_name} has {count} rows.")
    print("Total number of rows in the database: ", sum(table_counts.values()))


def calculateCoverage(session):
    # metabolites
    unique_cohort_ids_metabolite = session.query(distinct(CohortReferencesMetabolite.cohort_id)).count()
    unique_hmdb_ids_count = session.query(distinct(CohortMetabolite.cohort_id)).count()
    try:
        metabolite_coverage = round(unique_cohort_ids_metabolite / unique_hmdb_ids_count, 3)
    except:
        metabolite_coverage = "NA"
    # proteins

    unique_cohort_ids_protein = session.query(distinct(CohortReferencesProtein.cohort_id)).count()
    unique_uniprot_ids_count = session.query(distinct(CohortProtein.cohort_id)).count()
    try:
        protein_coverage = round(unique_cohort_ids_protein / unique_uniprot_ids_count, 3)
    except:
        protein_coverage = "NA"
    # phenotypes
    unique_cohort_ids_phenotype_hpo = session.query(distinct(CohortReferencesPhenotype.cohort_id)).count()
    unique_cohort_ids_phenotype_mondo = session.query(distinct(CohortReferencesDisease.cohort_id)).count()

    unique_mondo_ids_count = session.query(distinct(CohortReferencesDisease.cohort_id)).count()
    unique_hpo_ids_count = session.query(distinct(CohortPhenotype.cohort_id)).count()
    unique_hpo_mondo_ids_count = unique_cohort_ids_phenotype_hpo + unique_cohort_ids_phenotype_mondo
    try:
        phenotype_coverage_hpo = round(unique_cohort_ids_phenotype_hpo / unique_hpo_ids_count, 3)
        phenotype_coverage_mondo = round(unique_cohort_ids_phenotype_mondo / unique_hpo_ids_count, 3)
        phenotype_coverage_mondo_hpo = round(unique_hpo_mondo_ids_count / unique_hpo_mondo_ids_count,3)
    except:
        phenotype_coverage = "NA"
    # genomic_variant
    unique_cohort_ids_genomic_variant = session.query(distinct(CohortReferencesVariant.cohort_id)).count()
    unique_clinvar_ids_count = session.query(distinct(CohortGenomicVariant.cohort_id)).count()
    try:
        genomic_variant_coverage = round(unique_cohort_ids_genomic_variant / unique_clinvar_ids_count, 3)
    except:
        genomic_variant_coverage = "NA"

    # Output the results
    coverage = [
        ("Unique metabolite_cohort count", unique_cohort_ids_metabolite),
        ("Unique hmdb_id count", unique_hmdb_ids_count),
        ("metabolite_coverage", metabolite_coverage),

        ("Unique protein_cohort count", unique_cohort_ids_protein),
        ("Unique uniprot_id count", unique_uniprot_ids_count),
        ("protein_coverage", protein_coverage),

        ("Unique phenotype_cohort_hpo count", unique_cohort_ids_phenotype_hpo),
        ("Unique hpo_id count", unique_hpo_ids_count),
        ("phenotype_coverage_hpo", phenotype_coverage_hpo),
        ("Unique mondo_id count", unique_mondo_ids_count),
        ("phenotype_coverage_mondo", phenotype_coverage_mondo),
        ("phenotype_coverage_total", phenotype_coverage_mondo_hpo),

        ("Unique genomic variants_cohort_id count", unique_cohort_ids_genomic_variant),
        ("Unique clinvar_id count ", unique_clinvar_ids_count),
        ("genomic_variants", genomic_variant_coverage)
    ]
    coveragefilename = '../data/DyHealthNet/coverage_summary.csv'
    if not os.path.exists(os.path.dirname(coveragefilename)):
        os.makedirs(os.path.dirname(coveragefilename))
    with open(coveragefilename, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Writing the header
        writer.writerow(["ValueName", "Value"])
        # Writing the data
        writer.writerows(coverage)

    print(f"Data successfully written to {coveragefilename}")
    print(coverage)


def testingSetup(session):
    protein = test_protein(5)
    gene = test_gene()
    genomicVariant = test_variant()
    interactions = test_proteinAssocProtein(5)
    test = 2


if __name__ == '__main__':
    # Variant_affects_gene.__table__.drop(engine, checkfirst=True)
    # Define a session
    Session = sessionmaker(bind=engine)
    db_session = Session()
    #delete_tables(db_session)
    dotenv.load_dotenv()
    metadata = MetaData()
    create_tables()
    protein_data_path = PROTEIN_PATH
    pheno_data_path = PHENO_PATH
    metabo_data_path = METABOLITE_PATH
    edges_path = EDGES_PATH
    genomic_variant_meta_path = GENOMIC_VARIANT_META_PATH
    gwas_stats_path = GWAS_STATS_PATH
    data_dir = DATA_DIR

    if not all([pheno_data_path, protein_data_path, metabo_data_path,
                genomic_variant_meta_path, gwas_stats_path, edges_path]):
        raise ValueError("Please provide paths to the phenotype, protein, metabolite and edges files.")

    if not all([os.path.exists(x) for x in [pheno_data_path, protein_data_path, metabo_data_path,
                                            edges_path, genomic_variant_meta_path, gwas_stats_path]]):
        raise ValueError("Some of the provided paths do not exist.")

    add_genomic_variants_neddrex(db_session, genomic_variant_meta_path, obs_source=OBSERVATIONS)
    """
    add_disorder_data(db_session, pheno_data_path, obs_source=OBSERVATIONS)
    add_phenotype_data(db_session, pheno_data_path, obs_source=OBSERVATIONS, data_dir=data_dir)
    add_protein_data(db_session, protein_data_path, obs_source=OBSERVATIONS)
    add_metabolite_data(db_session, metabo_data_path, obs_source=OBSERVATIONS, data_dir=data_dir)

    # second pass for phenotypes
    add_phenotype_data(db_session, pheno_data_path, obs_source='external', data_dir=data_dir)

    # add cohort phenotype data as the mapping is incomplete
    
    add_cohort_phenotype_data(db_session, pheno_data_path, obs_source=OBSERVATIONS)
    add_cohort_metabolite_data(db_session, metabo_data_path, obs_source=OBSERVATIONS)
    add_cohort_protein_data(db_session, protein_data_path, obs_source=OBSERVATIONS)
    """
    add_cohort_genomic_variants(db_session, genomic_variant_meta_path, gwas_stats_path, obs_source=OBSERVATIONS)

    # add the edges calculated from the available data
    # add_calculated_edges(db_session, edges_path, pheno_data_path, protein_data_path, metabo_data_path)

    # calculateCoverage(db_session)
    # count the number of entries in the database
    metadata.reflect(bind=engine)
    # countEntries(db_session, metadata)

    # add remaining things (indexes, views)
    add_views(db_session)
    add_indexes(db_session, engine, metadata)
    db_session.close()
    print("Database setup complete.")
