from sqlalchemy import create_engine, Column, String, Integer, ForeignKey
from sqlalchemy.orm import declarative_base


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

