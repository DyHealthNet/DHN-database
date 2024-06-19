from sqlalchemy import Column, String, Integer, ForeignKey, ARRAY
from sqlalchemy.orm import declarative_base

# Create a declarative base
Base = declarative_base()


class Gene(Base):
    __tablename__ = 'genes'
    entrez_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    synonyms = Column(ARRAY(String))
    chromosome = Column(String)
    observation_source = Column(String)


class Disorder(Base):
    __tablename__ = 'disorders'
    mondo_id = Column(String, primary_key=True)
    description = Column(String)
    xrefs = Column(ARRAY(String))  # take from nedrex
    observation_source = Column(String)


class Phenotype(Base):
    __tablename__ = 'phenotypes'
    hpo_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    xrefs = Column(ARRAY(String))
    synonyms = Column(String)
    observation_source = Column(String)


class Protein(Base):
    __tablename__ = 'proteins'
    uniprot_id = Column(String, primary_key=True)
    sequence = Column(String)
    gene_entrez_id = Column(String)
    description = Column(String)
    observation_source = Column(String)


class Metabolite(Base):
    __tablename__ = 'metabolites'
    hmdb_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    synonyms = Column(String)
    xrefs = Column(ARRAY(String))
    observation_source = Column(String)


class GeneAssocDisorder(Base):
    __tablename__ = 'gene_associates_disorders'
    id = Column(Integer, primary_key=True)
    entrez_id = Column(String, ForeignKey('genes.entrez_id'))
    mondo_id = Column(String, ForeignKey('disorders.mondo_id'))
    edge_source = Column(String)


class DisorderAssocPhenotype(Base):
    __tablename__ = 'disorder_associates_phenotypes'
    id = Column(Integer, primary_key=True)
    mondo_id = Column(String, ForeignKey('disorders.mondo_id'))
    hpo_id = Column(String, ForeignKey('phenotypes.hpo_id'))
    edge_source = Column(String)


class ProteinAssocProtein(Base):
    __tablename__ = 'protein_associates_proteins'
    id = Column(Integer, primary_key=True, autoincrement=True)
    uniprot_id_memberOne = Column(String)
    uniprot_id_memberTwo = Column(String)


class ProteinAssocMetabolite(Base):
    __tablename__ = 'protein_associates_metabolites'
    id = Column(Integer, primary_key=True)
    uniprot_id = Column(String, ForeignKey('proteins.uniprot_id'))
    hmdb_id = Column(String, ForeignKey('metabolites.hmdb_id'))
    edge_source = Column(String)


class MetaboliteAssocDisorder(Base):
    __tablename__ = 'metabolite_associates_disorders'
    id = Column(Integer, primary_key=True)
    hmdb_id = Column(String, ForeignKey('metabolites.hmdb_id'))
    mondo_id = Column(String, ForeignKey('disorders.mondo_id'))
    edge_source = Column(String)
