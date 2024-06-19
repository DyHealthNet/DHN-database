import os

import networkx as nx
from sqlalchemy import URL, create_engine
from sqlalchemy.orm import sessionmaker, DeclarativeBase

from metabolite_mapping import read_metabolite_mapping, read_hmdb_data, download_metabolite_data, \
    retrieve_assoc_metabolite_nodes
from models import *
from query_nedrex import get_needed_snomed_ids, domain_id_to_mondo, get_disorder_data, get_edge_associations, \
    get_harmonizome_data, get_gene_data, get_phenotype_data
from hpo_mapping import download_hpo_ontology, read_hpo_ontology, ontology_data_to_network, snomed_from_hpo
from protein_mapping import get_proteinID_neddrex, read_proteinID_chris, add_proteinSet_data, \
    retrieve_interacting_proteins_neo4j, get_protein_nodes
from nedrex.core import get_collection_attributes
from sqlalchemy import create_engine, MetaData, create_engine, inspect, Table
from sqlalchemy.ext.declarative import declarative_base

# create postrgres db engine in memory
url = url_object = URL.create(
    "postgresql",
    username="postgres",
    password="password",  # plain (unescaped) text
    host="0.0.0.0",
    port=9000,
    database="postgres",
)
engine = create_engine(url)


def create_tables():
    # does not recreate tables if they already exist
    Base.metadata.create_all(engine)


def example_query(session):
    # Querying the database
    results = (session
               .query(Gene, Phenotype)
               .join(GeneAssocDisorder, Gene.entrez_id == GeneAssocDisorder.entrez_id)
               .filter(Gene.entrez_id == "12345")
               .all())
    for gene, phenotype in results:
        print(
            f"Gene entrez_id: {gene.entrez_id}, Phenotype hpo_id: {phenotype.hpo_id},"
            f" SNOMED ID: {phenotype.snomed_id}, OMIM ID: {phenotype.omim_id}")


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
                           xrefs: dict, gene_info: dict, assoc_graph: nx.Graph, obs_source: str = None) -> tuple[set, set, set, int]:
    """
    Queries the needed snomed ids and retrieves the associated genes and disorders from NeDRex
    :param gene_info: Information about the genes needed for the database (display name, synonyms, etc.)
    :param xrefs: cross references for the mondo ids to other databases
    :param descriptions: descriptions for the mondo ids
    :param needed_snomed: snomed ids in the dataset
    :param snomed_to_mondo: map from snomed to mondo ids from nedrex
    :param assoc_graph: association graph from nedrex of mondo ids to genes
    :param obs_source: Describes the source of observations - e.g. CHRIS
    :return: list of genes to add, list of disorders to add, list of gene associations to add, number of snomed ids found
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

            disorders.add(Disorder(mondo_id=mondo_id, xrefs=xref, description=description, observation_source=obs_source))
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


def add_items(session, items: iter, column: type[DeclarativeBase], filter_args: list):
    """
    Adds items to the database if they do not already exist
    :param session: Session object
    :param items: iterable of items to add
    :param column: the type of item to add, must be a class from models.py
    :param filter_args: the attributes to filter by to check if the item already exists
    :return: None
    """
    for item in items:
        filter_values = {key: getattr(item, key) for key in filter_args}
        exists = session.query(column).filter_by(**filter_values).first()
        if exists is not None:
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
    assoc_graph = get_edge_associations(set(domain_to_mondo.values()), edge_type='gene_associated_with_disorder')
    # assoc_graph is filtered for ids that we need, now we can get the data for all genes in the graph since they're
    # all associated with the mondo ids
    gene_info = get_gene_data()
    gene_dict = {x['primaryDomainId']: x for x in gene_info}

    mondo_description = {mondo: data[mondo]['description'] for mondo in domain_to_mondo.values()}

    genes_to_add, disorders, gene_associations, found = retrieve_disorder_data(needed_snomed, domain_to_mondo,
                                                                               mondo_description, xrefs, gene_dict,
                                                                               assoc_graph, obs_source)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    add_items(session, disorders, Disorder, ['mondo_id'])
    add_items(session, gene_associations, GeneAssocDisorder, ['entrez_id', 'mondo_id'])
    session.commit()
    print(f"Found and successfully added {found} snomed ids with diseases to db")


def add_phenotype_data(session, phenotype_path: str, data_dir: str = '../data', obs_source: str = None):
    """
    Adds phenotype data to the database given a path to a file with phenotype data
    :param obs_source: Describes the source of observations - e.g. CHRIS
    :param session: Database session object
    :param phenotype_path: str, path to file with phenotype data
    :param data_dir: str, path to the data directory
    :return: None
    """
    # data handling
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    needed_files = [f'{data_dir}/hp.json', f'{data_dir}/phenotype.hpoa']
    if not all([os.path.exists(f) for f in needed_files]):
        download_hpo_ontology(data_dir)
    needed_ids = get_needed_snomed_ids(phenotype_path)

    hpo_data = read_hpo_ontology(needed_files[0])
    hpo_graph = ontology_data_to_network(hpo_data)

    available_snomed_ids = snomed_from_hpo(hpo_graph, needed_ids)

    #         hpo_id = hpo_id.replace(':', '.').replace('HP', 'hpo')
    available_snomed_ids = {k: v.replace(':', '.').replace('HP', 'hpo') for k, v in available_snomed_ids.items()}
    pheno_data = get_phenotype_data(set(available_snomed_ids.values()))

    additional_data = {item['primaryDomainId']: item for item in pheno_data}
    genes_to_add, phenotypes, disorder_associations, _ = retrieve_phenotype_data(available_snomed_ids, additional_data, obs_source)

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


def add_protein_data(session, proteinData_path, observation_source):
    proteinIds = read_proteinID_chris(proteinData_path)
    proteinNodes = get_protein_nodes(proteinIds, observation_source)
    test = 2
    add_items(session, proteinNodes, Protein, ['uniprot_id'])
    session.commit()


def add_metabolite_data(session, metabolite_path, data_dir: str = '../data', obs_source: str = None):
    hmdb_data_path = f'{data_dir}/hmdb_metabolites.xml'
    download_metabolite_data(data_dir)
    metabolite_mapping = read_metabolite_mapping(metabolite_path)
    unique_metabolites = set()

    # split metabolites that have ; in them
    for metabolite in metabolite_mapping['hmdb_id'].dropna():
        unique_metabolites.update(metabolite.split(';'))

    print(f"Found {len(unique_metabolites)} unique metabolites in the mapping file.")

    hmdb_mapping = read_hmdb_data(hmdb_data_path, unique_metabolites)
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
    # add_missing(session, missing_proteins, 'proteins')

    for metabolite in hmdb_mapping:
        metabolite_name = f"hmdb.{metabolite}"

        metabolites.append(Metabolite(hmdb_id=metabolite_name, display_name=hmdb_mapping[metabolite]['display_name'],
                                      description=hmdb_mapping[metabolite]['description'],
                                      synonyms=hmdb_mapping[metabolite]['synonyms'],
                                      xrefs=hmdb_mapping[metabolite]['xrefs'],
                                      observation_source=obs_source))

        for disease in hmdb_mapping[metabolite]['diseases']:
            if session.query(Disorder).filter_by(mondo_id=omim_mapping.get(f"omim.{disease}", None)).first() is None:
                continue
            metabolite_disease_associations.append(MetaboliteAssocDisorder(hmdb_id=metabolite_name,
                                                                           mondo_id=omim_mapping[f"omim.{disease}"]))

        for protein in hmdb_mapping[metabolite]['proteins']:
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

    print(f"Got {len(data)} missing {node_type} to add to the database.")
    valid_node_types[node_type](session, missing_ids=data, obs_source='external')


if __name__ == '__main__':
    # cohort study
    observations = "CHRIS"
    # Define a session
    Session = sessionmaker(bind=engine)
    session = Session()

    create_tables()
    pheno_data_path = '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv'
    protein_data_path = '../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt'
    metabo_data_path = '../data/DyHealthNet/chris_summary_data/metabolites/CHRIS_biocristes7500SumStats.txt'
    inspector = inspect(engine)
    #metadata = MetaData()
    #metadata.reflect(bind=engine)
    #metadata.drop_all(bind=engine)
    inspector = inspect(engine)
    # Get and print the names of all tables
   # table_names = inspector.get_table_names()
   # print("Tables and their columns in the database:")
    #for table_name in table_names:
     #   print(f"\nTable: {table_name}")
        # Get the columns for each table
      #  columns = inspector.get_columns(table_name)
       # for column in columns:
        #    column_name = column['name']
         #   column_type = column['type']
          #  print(f" - Column: {column_name}, Type: {column_type}")
    #add_disorder_data(session, pheno_data_path, obs_source=observations)
    # add_phenotype_data(session, pheno_data_path, obs_source=observations)
    # add_protein_data(session, protein_data_path, obs_source=observations)
   # proteinIds = set(read_proteinID_chris(protein_data_path))
   # protein_test_ids = ['P51814', 'P19419', 'P43080', 'Q14457', 'Q01968', 'O95238']
    #print("test:" ,test)
    metadata = MetaData()
    #add_protein_data(session, protein_data_path,observation_source='CHRIS')
    specific_table = Table('proteins', metadata, autoload_with=engine)

    # Query the table
    query = session.query(specific_table).count()  # Limit to 10 entries for display
    test =2

    #test = 2
    #proteinNodes = get_protein_nodes(proteinIds)
    #nodeExample = {'primaryDomainId': 'uniprot.P31946',
    # 'comments': "FUNCTION: Adapter protein implicated in the regulation of a large spectrum of both general and specialized signaling pathways. Binds to a large number of partners, usually by recognition of a phosphoserine or phosphothreonine motif. Binding generally results in the modulation of the activity of the binding partner. Negative regulator of osteogenesis. Blocks the nuclear translocation of the phosphorylated form (by AKT1) of SRPK2 and antagonizes its stimulatory effect on cyclin D1 expression resulting in blockage of neuronal apoptosis elicited by SRPK2. Negative regulator of signaling cascades that mediate activation of MAP kinases via AKAP13. {ECO:0000269|PubMed:17717073, ECO:0000269|PubMed:19592491, ECO:0000269|PubMed:21224381}.\nSUBUNIT: Homodimer (PubMed:17717073). Interacts with SAMSN1 and PRKCE (By similarity). Interacts with AKAP13 (PubMed:21224381). Interacts with SSH1 and TORC2/CRTC2 (PubMed:15454081, PubMed:15159416). Interacts with ABL1; the interaction results in cytoplasmic location of ABL1 and inhibition of cABL-mediated apoptosis (PubMed:15696159). Interacts with ROR2 (dimer); the interaction results in phosphorylation of YWHAB on tyrosine residues (PubMed:17717073). Interacts with GAB2 (PubMed:19172738). Interacts with YAP1 (phosphorylated form) (PubMed:17974916). Interacts with the phosphorylated (by AKT1) form of SRPK2 (PubMed:19592491). Interacts with PKA-phosphorylated AANAT (PubMed:11427721). Interacts with MYO1C (PubMed:24636949). Interacts with SIRT2 (PubMed:18249187). Interacts with the 'Thr-369' phosphorylated form of DAPK2 (PubMed:26047703). Interacts with PI4KB, TBC1D22A and TBC1D22B (PubMed:23572552). Interacts with the 'Ser-1134' and 'Ser-1161' phosphorylated form of SOS1 (PubMed:22827337). Interacts (via phosphorylated form) with YWHAB; this interaction occurs in a protein kinase AKT1-dependent manner (PubMed:15538381). Interacts with SLITRK1 (PubMed:19640509). Interacts with SYNPO2 (phosphorylated form); YWHAB competes with ACTN2 for interaction with SYNPO2 (By similarity). Interacts with RIPOR2 (via phosphorylated form) isoform 2; this interaction occurs in a chemokine-dependent manner and does not compete for binding of RIPOR2 with RHOA nor blocks inhibition of RIPOR2- mediated RHOA activity (PubMed:25588844). Interacts with MARK2 and MARK3 (PubMed:16959763). Interacts with TESK1; the interaction is dependent on the phosphorylation of TESK1 'Ser-437' and inhibits TESK1 kinase activity (PubMed:11555644). Interacts with MEFV (PubMed:27030597). Interacts with HDAC4 (PubMed:33537682). Interacts with ADAM22 (via C-terminus) (PubMed:15882968). {ECO:0000250|UniProtKB:Q9CQV8, ECO:0000269|PubMed:11427721, ECO:0000269|PubMed:11555644, ECO:0000269|PubMed:15159416, ECO:0000269|PubMed:15454081, ECO:0000269|PubMed:15538381, ECO:0000269|PubMed:15696159, ECO:0000269|PubMed:15882968, ECO:0000269|PubMed:16959763, ECO:0000269|PubMed:17085597, ECO:0000269|PubMed:17717073, ECO:0000269|PubMed:17974916, ECO:0000269|PubMed:18249187, ECO:0000269|PubMed:19172738, ECO:0000269|PubMed:19592491, ECO:0000269|PubMed:19640509, ECO:0000269|PubMed:21224381, ECO:0000269|PubMed:22827337, ECO:0000269|PubMed:23572552, ECO:0000269|PubMed:24636949, ECO:0000269|PubMed:25588844, ECO:0000269|PubMed:26047703, ECO:0000269|PubMed:27030597, ECO:0000269|PubMed:33537682}.\nSUBUNIT: (Microbial infection) Interacts with herpes simplex virus 1 protein UL46. {ECO:0000269|PubMed:23938468}.\nSUBUNIT: (Microbial infection) Probably interacts with Chlamydia trachomatis protein IncG. {ECO:0000305|PubMed:11260479}.\nINTERACTION: P31946; Q9P0K1-3: ADAM22; NbExp=2; IntAct=EBI-359815, EBI-1567267; P31946; Q12802: AKAP13; NbExp=3; IntAct=EBI-359815, EBI-1373806; P31946; Q96B36: AKT1S1; NbExp=3; IntAct=EBI-359815, EBI-720593; P31946; A0A0S2Z5Q7: ALS2; NbExp=3; IntAct=EBI-359815, EBI-25928834; P31946; P05067: APP; NbExp=3; IntAct=EBI-359815, EBI-77613; P31946; P54253: ATXN1; NbExp=5; IntAct=EBI-359815, EBI-930964; P31946; Q92934: BAD; NbExp=5; IntAct=EBI-359815, EBI-700771; P31946; P15056: BRAF; NbExp=8; IntAct=EBI-359815, EBI-365980; P31946; P22681: CBL; NbExp=4; IntAct=EBI-359815, EBI-518228; P31946; O00257: CBX4; NbExp=3; IntAct=EBI-359815, EBI-722425; P31946; P30304: CDC25A; NbExp=10; IntAct=EBI-359815, EBI-747671; P31946; P30305: CDC25B; NbExp=5; IntAct=EBI-359815, EBI-1051746; P31946; P30307: CDC25C; NbExp=6; IntAct=EBI-359815, EBI-974439; P31946; O94921: CDK14; NbExp=6; IntAct=EBI-359815, EBI-1043945; P31946; Q53ET0: CRTC2; NbExp=4; IntAct=EBI-359815, EBI-1181987; P31946; Q9NYF0: DACT1; NbExp=4; IntAct=EBI-359815, EBI-3951744; P31946; Q13627-2: DYRK1A; NbExp=3; IntAct=EBI-359815, EBI-1053621; P31946; Q9UQC2: GAB2; NbExp=6; IntAct=EBI-359815, EBI-975200; P31946; P55040: GEM; NbExp=4; IntAct=EBI-359815, EBI-744104; P31946; P56524: HDAC4; NbExp=5; IntAct=EBI-359815, EBI-308629; P31946; P42858: HTT; NbExp=11; IntAct=EBI-359815, EBI-466029; P31946; Q5S007: LRRK2; NbExp=5; IntAct=EBI-359815, EBI-5323863; P31946; Q99759: MAP3K3; NbExp=4; IntAct=EBI-359815, EBI-307281; P31946; Q99683: MAP3K5; NbExp=3; IntAct=EBI-359815, EBI-476263; P31946; Q7KZI7: MARK2; NbExp=5; IntAct=EBI-359815, EBI-516560; P31946; P27448: MARK3; NbExp=8; IntAct=EBI-359815, EBI-707595; P31946; P26045: PTPN3; NbExp=6; IntAct=EBI-359815, EBI-1047946; P31946; P04049: RAF1; NbExp=35; IntAct=EBI-359815, EBI-365996; P31946; Q96TC7: RMDN3; NbExp=5; IntAct=EBI-359815, EBI-1056589; P31946; P61587: RND3; NbExp=2; IntAct=EBI-359815, EBI-1111534; P31946; Q96JI7: SPG11; NbExp=2; IntAct=EBI-359815, EBI-2822128; P31946; P78362: SRPK2; NbExp=2; IntAct=EBI-359815, EBI-593303; P31946; Q8WYL5: SSH1; NbExp=3; IntAct=EBI-359815, EBI-1222387; P31946; P49815: TSC2; NbExp=6; IntAct=EBI-359815, EBI-396587; P31946; P46937: YAP1; NbExp=7; IntAct=EBI-359815, EBI-1044059; P31946; P31946: YWHAB; NbExp=3; IntAct=EBI-359815, EBI-359815; P31946; P62258: YWHAE; NbExp=12; IntAct=EBI-359815, EBI-356498; P31946; P61981: YWHAG; NbExp=5; IntAct=EBI-359815, EBI-359832; P31946; P27348: YWHAQ; NbExp=5; IntAct=EBI-359815, EBI-359854; P31946; P67828: CSNK1A1; Xeno; NbExp=3; IntAct=EBI-359815, EBI-7540603; P31946; P55041: Gem; Xeno; NbExp=3; IntAct=EBI-359815, EBI-7082069; P31946; Q11184: let-756; Xeno; NbExp=2; IntAct=EBI-359815, EBI-3843983; P31946; P61588: Rnd3; Xeno; NbExp=5; IntAct=EBI-359815, EBI-6930266; P31946; Q91YE8: Synpo2; Xeno; NbExp=3; IntAct=EBI-359815, EBI-7623057; P31946; B7UM99: tir; Xeno; NbExp=2; IntAct=EBI-359815, EBI-2504426; P31946; P22893: Zfp36; Xeno; NbExp=5; IntAct=EBI-359815, EBI-647803; P31946; Q76353; Xeno; NbExp=3; IntAct=EBI-359815, EBI-6248077;\nSUBCELLULAR LOCATION: Cytoplasm {ECO:0000269|PubMed:17081065}. Melanosome {ECO:0000269|PubMed:17081065}. Note=Identified by mass spectrometry in melanosome fractions from stage I to stage IV.\nSUBCELLULAR LOCATION: Vacuole membrane {ECO:0000269|PubMed:11260479}. Note=(Microbial infection) Upon infection with Chlamydia trachomatis, this protein is associated with the pathogen-containing vacuole membrane where it colocalizes with IncG. {ECO:0000269|PubMed:11260479}.\nALTERNATIVE PRODUCTS: Event=Alternative initiation; Named isoforms=2; Name=Long;   IsoId=P31946-1; Sequence=Displayed; Name=Short;   IsoId=P31946-2; Sequence=VSP_018632;\nPTM: The alpha, brain-specific form differs from the beta form in being phosphorylated. Phosphorylated on Ser-60 by protein kinase C delta type catalytic subunit in a sphingosine-dependent fashion. {ECO:0000250}.\nSIMILARITY: Belongs to the 14-3-3 family. {ECO:0000305}.",
    # 'created': '2024-04-08T13:50:46.372000', 'dataSources': ['uniprot'], 'displayName': '1433B_HUMAN',
    # 'domainIds': ['uniprot.P31946', 'ensembl.ENSP00000300161.4', 'ensembl.ENSP00000361930.3'], 'geneName': 'YWHAB',
     #'sequence': 'MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFYLKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGDAGEGEN',
     #'synonyms': ['14-3-3 protein beta/alpha', 'Protein 1054', 'Protein kinase C inhibitor protein 1', 'KCIP-1'],
     #'taxid': 9606, 'type': 'Protein', 'updated': '2024-04-08T13:50:46.372000'}
    #print(nodeExample.get('primaryDomainId'))
    #print(nodeExample['domainIds'][0])

    #testProt = Protein(uniprot_id=nodeExample['domainIds'[1]])


'''  class Protein(Base):
        __tablename__ = 'proteins'
        uniprot_id = Column(String, primary_key=True)
        sequence = Column(String)
        gene_entrez_id = Column(String, ForeignKey('genes.entrez_id'))
        description = Column(String)
        observation_source = Column(String)
'''
    # add_metabolite_data(session, metabo_data_path, obs_source=observations)
    #retrieve_interacting_proteins_neo4j(protein_test_ids)
