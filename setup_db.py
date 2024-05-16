from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from models import Base, Gene, Disorder, Phenotype, GeneAssocDisorder, GeneAssocPhenotype

# create postrgres db engine in memory
engine = create_engine("postgresql://postgres:password@172.17.0.2:5432/postgres")


def example_create():
    # does not recreate tables if they already exist
    Base.metadata.create_all(engine)


def example_add(session):
    # example test
    new_gene = Gene(entrezid=12345)
    new_disorder = Disorder(mondoid='MONDO:0000001', snomed_id='SNOMED:0000001')
    new_phenotype = Phenotype(hpoid='HP:0000001', snomed_id='SNOMED:0000002', omim_id='OMIM:0000001')

    session.add(new_gene)
    session.add(new_disorder)
    session.add(new_phenotype)
    session.commit()

    # Create associations
    gene_assoc_disorder = GeneAssocDisorder(entrezid=12345, mondoid='MONDO:0000001', edge_source='source1')
    gene_assoc_phenotype = GeneAssocPhenotype(entrezid=12345, hpoid='HP:0000001')

    session.add(gene_assoc_disorder)
    session.add(gene_assoc_phenotype)
    session.commit()


def example_query(session):
    # Querying the database
    results = (session
               .query(Gene, Phenotype)
               .join(GeneAssocPhenotype, Gene.entrezid == GeneAssocPhenotype.entrezid)
               .join(Phenotype, GeneAssocPhenotype.hpoid == Phenotype.hpoid)
               .filter(Gene.entrezid == 12345)
               .all())
    for gene, phenotype in results:
        print(
            f"Gene EntrezID: {gene.entrezid}, Phenotype HPOID: {phenotype.hpoid},"
            f" SNOMED ID: {phenotype.snomed_id}, OMIM ID: {phenotype.omim_id}")


if __name__ == '__main__':
    # Define a session
    Session = sessionmaker(bind=engine)
    session = Session()

    example_create()
    example_add(session)
    example_query(session)
