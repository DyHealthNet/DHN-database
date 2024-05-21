from sqlalchemy import create_engine, Column, String, Integer, ForeignKey
from sqlalchemy.orm import declarative_base


# Create a declarative base
Base = declarative_base()


class Gene(Base):
    __tablename__ = 'genes'
    entrez_id = Column(String, primary_key=True)


class Disorder(Base):
    __tablename__ = 'disorders'
    mondo_id = Column(String, primary_key=True)
    snomed_id = Column(String)


class Phenotype(Base):
    __tablename__ = 'phenotypes'
    hpo_id = Column(String, primary_key=True)
    snomed_id = Column(String, nullable=True)
    omim_id = Column(String)
    orpha_id = Column(String)


class GeneAssocDisorder(Base):
    __tablename__ = 'gene_associates_disorders'
    id = Column(Integer, primary_key=True)
    entrez_id = Column(String, ForeignKey('genes.entrez_id'))
    mondo_id = Column(String, ForeignKey('disorders.mondo_id'))
    edge_source = Column(String)


class GeneAssocPhenotype(Base):
    __tablename__ = 'gene_associates_phenotypes'
    id = Column(Integer, primary_key=True)
    entrez_id = Column(String, ForeignKey('genes.entrez_id'))
    hpo_id = Column(String, ForeignKey('phenotypes.hpo_id'))
    edge_source = Column(String)


class Protein(Base):
    __tablename__ = 'proteins'
    id = Column(Integer, primary_key=True)
    uniprotId = Column(String, primary_key=True)
    gene = Column(String, ForeignKey('genes.entrez_id'))


