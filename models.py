from sqlalchemy import Column, String, Integer, ForeignKey, ARRAY, Float
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


# Association tables between node types
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


# Effects tables based on the calculations done on the cohort study
class EffectsProteinProtein(Base):
    __tablename__ = 'effects_protein_protein'
    id = Column(Integer, primary_key=True, autoincrement=True)
    uniprot_id_1 = Column(String, ForeignKey('proteins.uniprot_id'))
    uniprot_id_2 = Column(String, ForeignKey('proteins.uniprot_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsProteinMetabolite(Base):
    __tablename__ = 'effects_protein_metabolite'
    id = Column(Integer, primary_key=True, autoincrement=True)
    uniprot_id = Column(String, ForeignKey('proteins.uniprot_id'))
    hmdb_id = Column(String, ForeignKey('metabolites.hmdb_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsProteinPhenotype(Base):
    __tablename__ = 'effects_protein_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    uniprot_id = Column(String, ForeignKey('proteins.uniprot_id'))
    hpo_id = Column(String, ForeignKey('phenotypes.hpo_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsProteinDisorder(Base):
    __tablename__ = 'effects_protein_disorder'
    id = Column(Integer, primary_key=True, autoincrement=True)
    uniprot_id = Column(String, ForeignKey('proteins.uniprot_id'))
    mondo_id = Column(String, ForeignKey('disorders.mondo_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsMetaboliteMetabolite(Base):
    __tablename__ = 'effects_metabolite_metabolite'
    id = Column(Integer, primary_key=True, autoincrement=True)
    hmdb_id_1 = Column(String, ForeignKey('metabolites.hmdb_id'))
    hmdb_id_2 = Column(String, ForeignKey('metabolites.hmdb_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsMetabolitePhenotype(Base):
    __tablename__ = 'effects_metabolite_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    hmdb_id = Column(String, ForeignKey('metabolites.hmdb_id'))
    hpo_id = Column(String, ForeignKey('phenotypes.hpo_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsMetaboliteDisorder(Base):
    __tablename__ = 'effects_metabolite_disorder'
    id = Column(Integer, primary_key=True, autoincrement=True)
    hmdb_id = Column(String, ForeignKey('metabolites.hmdb_id'))
    mondo_id = Column(String, ForeignKey('disorders.mondo_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsPhenotypePhenotype(Base):
    __tablename__ = 'effects_phenotype_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    hpo_id_1 = Column(String, ForeignKey('phenotypes.hpo_id'))
    hpo_id_2 = Column(String, ForeignKey('phenotypes.hpo_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsPhenotypeDisorder(Base):
    __tablename__ = 'effects_phenotype_disorder'
    id = Column(Integer, primary_key=True, autoincrement=True)
    hpo_id = Column(String, ForeignKey('phenotypes.hpo_id'))
    mondo_id = Column(String, ForeignKey('disorders.mondo_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsDisorderDisorder(Base):
    __tablename__ = 'effects_disorder_disorder'
    id = Column(Integer, primary_key=True, autoincrement=True)
    mondo_id_1 = Column(String, ForeignKey('disorders.mondo_id'))
    mondo_id_2 = Column(String, ForeignKey('disorders.mondo_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)
