from sqlalchemy import create_engine, Column, String, Integer, ForeignKey
from sqlalchemy.orm import sessionmaker, relationship, declarative_base

# create postrgres db engine in memory
engine = create_engine("postgresql://postgres:password@172.17.0.2:5432/postgres")

# Create a declarative base
Base = declarative_base()


class Gene(Base):
    __tablename__ = 'genes'
    entrezid = Column(Integer, primary_key=True)


class Disorder(Base):
    __tablename__ = 'disorders'
    mondoid = Column(String, primary_key=True)
    snomed_id = Column(String)


class Phenotype(Base):
    __tablename__ = 'phenotypes'
    hpoid = Column(String, primary_key=True)
    snomed_id = Column(String)
    omim_id = Column(String)


class GeneAssocDisorder(Base):
    __tablename__ = 'gene_assoc_disorders'
    id = Column(Integer, primary_key=True)
    entrezid = Column(Integer, ForeignKey('genes.entrezid'))
    mondoid = Column(String, ForeignKey('disorders.mondoid'))
    edge_source = Column(String)


class GeneAssocPhenotype(Base):
    __tablename__ = 'gene_assoc_phenotypes'
    id = Column(Integer, primary_key=True)
    entrezid = Column(Integer, ForeignKey('genes.entrezid'))
    hpoid = Column(String, ForeignKey('phenotypes.hpoid'))


# does not recreate tables if they already exist
Base.metadata.create_all(engine)


def example_add():
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


if __name__ == '__main__':
    # Define a session
    Session = sessionmaker(bind=engine)
    session = Session()

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