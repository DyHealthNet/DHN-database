import os
from sqlalchemy import create_engine, URL
from sqlalchemy.orm import sessionmaker
from models import *
from query_nedrex import needed_snomed_ids, domain_id_to_mondo, get_disorder_data, get_all_associations, \
    get_harmonizome_data
from hpo_mapping import download_hpo_ontology, read_hpo_ontology, ontology_data_to_network, snomed_from_hpo, \
    hpo_to_xref, disorder_to_mondo
from protein_mapping import get_proteinID_neddrex, read_proteinID_chris

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



def example_add(session):
    # example test
    new_gene = Gene(entrez_id="12345")
    new_disorder = Disorder(mondo_id='MONDO:0000001', snomed_id='SNOMED:0000001')
    new_phenotype = Phenotype(hpo_id='HP:0000001', snomed_id='SNOMED:0000002', omim_id='OMIM:0000001')

    session.add(new_gene)
    session.add(new_disorder)
    session.add(new_phenotype)
    session.commit()

    # Create associations
    gene_assoc_disorder = GeneAssocDisorder(entrez_id="12345", mondo_id='MONDO:0000001', edge_source='source1')
    gene_assoc_phenotype = GeneAssocPhenotype(entrez_id="12345", hpo_id='HP:0000001')

    session.add(gene_assoc_disorder)
    session.add(gene_assoc_phenotype)
    session.commit()


def example_query(session):
    # Querying the database
    results = (session
               .query(Gene, Phenotype)
               .join(GeneAssocPhenotype, Gene.entrez_id == GeneAssocPhenotype.entrez_id)
               .join(Phenotype, GeneAssocPhenotype.hpo_id == Phenotype.hpo_id)
               .filter(Gene.entrez_id == "12345")
               .all())
    for gene, phenotype in results:
        print(
            f"Gene entrez_id: {gene.entrez_id}, Phenotype hpo_id: {phenotype.hpo_id},"
            f" SNOMED ID: {phenotype.snomed_id}, OMIM ID: {phenotype.omim_id}")


def mondo_in_association_graph(mondo_id, assoc_graph):
    """
    Retrieves the genes associated with a mondo id from the association graph
    :param mondo_id: A mondo id, can be None
    :param assoc_graph: The association graph from NEDRex
    :return: genes associated with the mondo id, source of the data, harmonizome data if not in the association graph
    """
    harm = None
    if mondo_id is None:
        return None, None, None
    if mondo_id in assoc_graph:
        genes = list(assoc_graph[mondo_id])
        source = 'nedrex'
    else:
        harm = get_harmonizome_data(mondo_id)
        if harm is None:
            return None, None, None
        genes = [gene for genes in harm.values() for gene in genes]
        source = 'harmonizome'
    return genes, source, harm


def retrieve_disorder_data(needed_snomed, snomed_to_mondo, assoc_graph):
    """
    Queries the needed snomed ids and retrieves the associated genes and disorders from NEDRex
    :param needed_snomed: snomed ids in the dataset
    :param snomed_to_mondo: map from snomed to mondo ids from nedrex
    :param assoc_graph: association graph from nedrex of mondo ids to genes
    :return: list of genes to add, list of disorders to add, list of gene associations to add, number of snomed ids found
    """
    gene_associations = set()
    genes_to_add = set()
    disorders = set()
    found = 0
    for snomed_id in needed_snomed['snomed_id'].unique():
        # check if ; in snomed id and if so, do this for all ids
        snomed_ids = str(snomed_id).split(';')
        for snomed_id in snomed_ids:
            mondo_id = snomed_to_mondo.get(snomed_id)
            genes, source, harm = mondo_in_association_graph(mondo_id, assoc_graph)
            if genes is None:
                continue
            # add genes to set
            genes_to_add.update([Gene(entrez_id=x) for x in genes])
            disorders.add(Disorder(mondo_id=mondo_id, snomed_id=snomed_id))
            # add gene associations to set for each source
            if source == 'nedrex':
                gene_associations.update([GeneAssocDisorder(entrez_id=gene, mondo_id=mondo_id, edge_source='nedrex')
                                          for gene in genes])
            elif source == 'harmonizome':
                for source in harm:
                    gene_associations.update([GeneAssocDisorder(entrez_id=gene, mondo_id=mondo_id, edge_source=source)
                                              for gene in harm[source]])
            found += 1
        continue
    return genes_to_add, disorders, gene_associations, found


def retrieve_phenotype_data(needed_snomed, hpo_graph, data_dir):
    """
    Retrieve the phenotype data for the needed snomed ids
    # HPO conversion: Pathway
    # HPO data (HPO_ID ----> SNOMED_ID) - look for needed SNOMED IDs
    # -> map to external database (SNOMED_ID -- HPO_ID --> OMIM_ID/ORPHA_ID)
    # -> map to Mondo (SNOMED_ID -- OMIM_ID/ORPHA_ID --> Mondo_ID)
    # -> get associated genes (SNOMED_ID -- Mondo_ID --> Genes)

    :param needed_snomed: set of snomed ids that are needed
    :param hpo_graph: Graph of the HPO ontology
    :return: dictionary with the phenotype data
    """
    found = 0
    genes_to_add = set()
    phenotypes = set()
    gene_associations = set()
    needed_ids = set(needed_snomed['snomed_id'].unique())
    # go through all the nodes in the HPO graph and find the ones that have xrefs to SNOMED
    available_snomed_ids = snomed_from_hpo(hpo_graph, needed_ids)
    # convert the snomed ids to OMIM ids
    snomed_to_xref = hpo_to_xref(data_dir, available_snomed_ids)
    # get the disorder data
    disorder_data = get_disorder_data()
    final_mapping = disorder_to_mondo(disorder_data, snomed_to_xref)

    # find the genes that are associated with the mondo ids
    assoc_graph = get_all_associations()
    for snomed_id, mondo_id in final_mapping.items():
        genes, source, harm = mondo_in_association_graph(mondo_id, assoc_graph)
        if genes is None:
            continue

        genes_to_add.update([Gene(entrez_id=x) for x in genes])

        hpo_id = available_snomed_ids[snomed_id]
        xref_id = snomed_to_xref[snomed_id]
        if xref_id.startswith('OMIM'):
            phenotypes.add(Phenotype(hpo_id=hpo_id, snomed_id=snomed_id, omim_id=xref_id))
        else:
            phenotypes.add(Phenotype(hpo_id=hpo_id, snomed_id=snomed_id, orpha_id=xref_id))

        if source == 'nedrex':
            gene_associations.update([GeneAssocPhenotype(entrez_id=gene, hpo_id=hpo_id, edge_source='nedrex')
                                      for gene in genes])
        elif source == 'harmonizome':
            for source in harm:
                gene_associations.update([GeneAssocPhenotype(entrez_id=gene, hpo_id=hpo_id, edge_source=source)
                                          for gene in harm[source]])

        found += 1
    return genes_to_add, phenotypes, gene_associations, found


def add_items(session, items: iter, column: type[Gene | Phenotype | Disorder | GeneAssocDisorder | GeneAssocPhenotype | Protein],
              filter_args: list):
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
        except Exception as e:
            print("Exception: ", e)
    session.commit()


def add_disorder_data(session, snomed_id_path: str):
    """
    Adds disorder data to the database given a path to a file with snomed ids
    :param session: Database session object
    :param snomed_id_path: str, path to file with snomed ids
    :return: None
    """
    needed_snomed = needed_snomed_ids(snomed_id_path)
    data = get_disorder_data()
    snomed_to_mondo = domain_id_to_mondo(data)
    assoc_graph = get_all_associations()

    genes_to_add, disorders, gene_associations, found = retrieve_disorder_data(needed_snomed, snomed_to_mondo, assoc_graph)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    add_items(session, disorders, Disorder, ['mondo_id'])
    add_items(session, gene_associations, GeneAssocDisorder, ['entrez_id', 'mondo_id'])
    session.commit()
    print(f"Found and successfully added {found} snomed ids with diseases to db")


def add_phenotype_data(session, phenotype_path: str, data_dir: str = '../data'):
    """
    Adds phenotype data to the database given a path to a file with phenotype data
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
    needed_ids = needed_snomed_ids(phenotype_path)

    hpo_data = read_hpo_ontology(needed_files[0])
    hpo_graph = ontology_data_to_network(hpo_data)
    genes_to_add, disorders, gene_associations, found = retrieve_phenotype_data(needed_ids, hpo_graph, data_dir)

    add_items(session, genes_to_add, Gene, ['entrez_id'])
    add_items(session, disorders, Phenotype, ['hpo_id'])
    add_items(session, gene_associations, GeneAssocPhenotype, ['entrez_id', 'hpo_id'])
    session.commit()
    print(f"Found and successfully added {found} snomed ids with phenotypes to db")
def add_protein_data(session, proteinData_path):
    proteinData =  read_proteinID_chris(proteinData_path)
    proteinGeneDict = {}
    for proteinID in proteinData:
        get_proteinID_neddrex(proteinID)[0]['geneName']
        add_items(session, Protein, proteinID)
        proteinGeneDict[proteinID] = get_proteinID_neddrex(proteinID)[0]['geneName']
    if (proteinID != None):
        add_items(session,proteinGeneDict, Protein)


if __name__ == '__main__':
    # Define a session
    Session = sessionmaker(bind=engine)
    session = Session()
    create_tables()
    #paths
    pheno_data_path = '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv'
    protein_data_path = '../data/DyHealthNet/chris_summary_data/proteins/CHRIS_somalogic_descriptive_statistic.txt'
 #   add_disorder_data(session, pheno_data_path)
    add_phenotype_data(session, pheno_data_path)
    add_protein_data(session, protein_data_path)


