from sqlalchemy import Column, String, Integer, ForeignKey, ARRAY, Float, Index
from sqlalchemy.orm import declarative_base

# Create a declarative base
Base = declarative_base()

### Tables of data from the cohort study ###

class CohortPhenotype(Base):
    __tablename__ = 'cohort_phenotype'
    cohort_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    mondo_id = Column(String, ForeignKey('disorder.mondo_id'))
    hpo_id = Column(String, ForeignKey('phenotype.hpo_id'))


class CohortProtein(Base):
    __tablename__ = 'cohort_protein'
    cohort_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    uniprot_id = Column(String, ForeignKey('protein.uniprot_id'))


class CohortMetabolite(Base):
    __tablename__ = 'cohort_metabolite'
    cohort_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    hmdb_id = Column(String, ForeignKey('metabolite.hmdb_id'))


### Tables of data from the external knowledge graph ###

class Gene(Base):
    __tablename__ = 'gene'
    entrez_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    synonyms = Column(ARRAY(String))
    chromosome = Column(String)
    observation_source = Column(String)


class Disorder(Base):
    __tablename__ = 'disorder'
    mondo_id = Column(String, primary_key=True)
    description = Column(String)
    xrefs = Column(ARRAY(String))  # take from nedrex
    observation_source = Column(String)


class Phenotype(Base):
    __tablename__ = 'phenotype'
    hpo_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    xrefs = Column(ARRAY(String))
    synonyms = Column(String)
    observation_source = Column(String)


class Protein(Base):
    __tablename__ = 'protein'
    uniprot_id = Column(String, primary_key=True)
    sequence = Column(String)
    gene_entrez_id = Column(String)
    description = Column(String)
    observation_source = Column(String)


class Metabolite(Base):
    __tablename__ = 'metabolite'
    hmdb_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    synonyms = Column(String)
    xrefs = Column(ARRAY(String))
    observation_source = Column(String)


class Genomic_variant(Base):
    __tablename__ = "genomic_variant"
    variant_primaryDomainId = Column(String, primary_key=True)  #clinvar.17735
    alternativeSequence = Column(String)  #'T',
    chromosome = Column(String)  # 'NW_009646201.1',
    created = Column(String)  # '2024-06-17T12:36:21.275000'
    dataSources = Column(String)  #['clinvar'],
    domainIds = Column(String)  #['clinvar.17735', 'dbsnp.1556058284']
    position = Column(String)  #83614,
    referenceSequence = Column(String)  # 'TC',
    type = Column(String)  # 'GenomicVariant'
    variantType = Column(String)  #'Deletion'}


### Association tables ###

# Calculated effects of cohort observations

class EffectsProteinProtein(Base):
    __tablename__ = 'effects_protein_protein'
    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id_1 = Column(String, ForeignKey('cohort_protein.cohort_id'))
    protein_id_2 = Column(String, ForeignKey('cohort_protein.cohort_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsProteinMetabolite(Base):
    __tablename__ = 'effects_protein_metabolite'
    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(String, ForeignKey('cohort_protein.cohort_id'))
    metabolite_id = Column(String, ForeignKey('cohort_metabolite.cohort_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsProteinPhenotype(Base):
    __tablename__ = 'effects_protein_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(String, ForeignKey('cohort_protein.cohort_id'))
    phenotype_id = Column(String, ForeignKey('cohort_phenotype.cohort_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsMetaboliteMetabolite(Base):
    __tablename__ = 'effects_metabolite_metabolite'
    id = Column(Integer, primary_key=True, autoincrement=True)
    phenotype_id_1 = Column(String, ForeignKey('cohort_metabolite.cohort_id'))
    phenotype_id_2 = Column(String, ForeignKey('cohort_metabolite.cohort_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsMetabolitePhenotype(Base):
    __tablename__ = 'effects_metabolite_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    metabolite_id = Column(String, ForeignKey('cohort_metabolite.cohort_id'))
    phenotype_id = Column(String, ForeignKey('cohort_phenotype.cohort_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


class EffectsPhenotypePhenotype(Base):
    __tablename__ = 'effects_phenotype_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    phenotype_id_1 = Column(String, ForeignKey('cohort_phenotype.cohort_id'))
    phenotype_id_2 = Column(String, ForeignKey('cohort_phenotype.cohort_id'))
    p_value = Column(Float)
    adjusted_p_value = Column(Float)
    effect_size = Column(Float)
    effect_size_type = Column(String)


# Associations between entities in the external knowledge graph
class GeneAssocDisorder(Base):
    __tablename__ = 'gene_associates_disorder'
    id = Column(Integer, primary_key=True)
    entrez_id = Column(String, ForeignKey('gene.entrez_id'))
    mondo_id = Column(String, ForeignKey('disorder.mondo_id'))
    edge_source = Column(String)


class DisorderAssocPhenotype(Base):
    __tablename__ = 'disorder_associates_phenotype'
    id = Column(Integer, primary_key=True)
    mondo_id = Column(String, ForeignKey('disorder.mondo_id'))
    hpo_id = Column(String, ForeignKey('phenotype.hpo_id'))
    edge_source = Column(String)


class ProteinAssocProtein(Base):
    __tablename__ = 'protein_associates_protein'
    id = Column(Integer, primary_key=True, autoincrement=True)
    uniprot_id_1 = Column(String, ForeignKey('protein.uniprot_id'))
    uniprot_id_2 = Column(String, ForeignKey('protein.uniprot_id'))


class Variant_affects_gene(Base):
    __tablename__ = 'variant_affects_gene'
    id = Column(Integer, primary_key=True, autoincrement=True)
    genomic_variant = Column(String, ForeignKey('genomic_variant.variant_primaryDomainId'))
    entrez_id = Column(String, ForeignKey('gene.entrez_id'))


class ProteinAssocMetabolite(Base):
    __tablename__ = 'protein_associates_metabolite'
    id = Column(Integer, primary_key=True)
    uniprot_id = Column(String, ForeignKey('protein.uniprot_id'))
    hmdb_id = Column(String, ForeignKey('metabolite.hmdb_id'))
    edge_source = Column(String)


class MetaboliteAssocDisorder(Base):
    __tablename__ = 'metabolite_associates_disorder'
    id = Column(Integer, primary_key=True)
    hmdb_id = Column(String, ForeignKey('metabolite.hmdb_id'))
    mondo_id = Column(String, ForeignKey('disorder.mondo_id'))
    edge_source = Column(String)
