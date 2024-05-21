from sqlalchemy import create_engine, URL
from sqlalchemy.orm import sessionmaker
from models import *
from query_nedrex import needed_snomed_ids, domain_id_to_mondo, get_disorder_data, get_all_associations, \
    get_harmonizome_data

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
            if mondo_id is None:
                continue
            if mondo_id in assoc_graph:
                genes = list(assoc_graph[mondo_id])
                source = 'nedrex'
            else:
                harm = get_harmonizome_data(mondo_id)
                if harm is None:
                    continue
                genes = [gene for genes in harm.values() for gene in genes]
                source = 'harmonizome'
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


if __name__ == '__main__':
    # Define a session
    Session = sessionmaker(bind=engine)
    session = Session()
    create_tables()
    add_disorder_data(session, '../data/DyHealthNet/chris_summary_data/phenotypes/pheno_meta_all.tsv')
