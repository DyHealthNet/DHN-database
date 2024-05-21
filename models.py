from sqlalchemy import create_engine, Column, String, Integer, ForeignKey
from sqlalchemy.orm import declarative_base


# Create a declarative base
Base = declarative_base()


class Gene(Base):
    __tablename__ = 'genes'
    entrezid = Column(String, primary_key=True)


class Disorder(Base):
    __tablename__ = 'disorders'
    mondoid = Column(String, primary_key=True)
    snomed_id = Column(String)


class Phenotype(Base):
    __tablename__ = 'phenotypes'
    hpoid = Column(String, primary_key=True)
    snomed_id = Column(String, nullable=True)
    omim_id = Column(String)
    orpha_id = Column(String)


class GeneAssocDisorder(Base):
    __tablename__ = 'gene_assoc_disorders'
    id = Column(Integer, primary_key=True)
    entrezid = Column(String, ForeignKey('genes.entrezid'))
    mondoid = Column(String, ForeignKey('disorders.mondoid'))
    edge_source = Column(String)


class GeneAssocPhenotype(Base):
    __tablename__ = 'gene_assoc_phenotypes'
    id = Column(Integer, primary_key=True)
    entrezid = Column(String, ForeignKey('genes.entrezid'))
    hpoid = Column(String, ForeignKey('phenotypes.hpoid'))

class Protein(Base):
    __tablename__ = 'proteins'
    id = Column(Integer, primary_key=True)
    uniprotId = Column(String, primary_key=True)
    gene = Column(String , ForeignKey('genes.entrezid'))


