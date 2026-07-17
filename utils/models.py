from sqlalchemy import Column, String, Integer, ForeignKey, ARRAY, Float, Index, JSON, UniqueConstraint
from sqlalchemy.orm import declarative_base

# Create a declarative base
Base = declarative_base()


### Tables of data from the cohort study ###
class CohortVariant(Base):
    __tablename__ = 'cohort_variant'
    cohort_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    xrefs = Column(String)


class CohortPhenotype(Base):
    __tablename__ = 'cohort_phenotype'
    cohort_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    xrefs = Column(String)


class CohortProtein(Base):
    __tablename__ = 'cohort_protein'
    cohort_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    xrefs = Column(String)


class CohortMetabolite(Base):
    __tablename__ = 'cohort_metabolite'
    cohort_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    xrefs = Column(String)


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
    display_name = Column(String)
    description = Column(String)
    xrefs = Column(ARRAY(String))
    observation_source = Column(String)


class Phenotype(Base):
    __tablename__ = 'phenotype'
    hpo_id = Column(String, primary_key=True)
    display_name = Column(String)
    description = Column(String)
    xrefs = Column(ARRAY(String))
    synonyms = Column(String)
    observation_source = Column(String)


# Read Data from file into Object Datastructure.

class Protein(Base):
    __tablename__ = 'protein'
    uniprot_id = Column(String, primary_key=True)
    display_name = Column(String)
    sequence = Column(String)
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


class GenomicVariant(Base):
    __tablename__ = "genomic_variant"
    clinvar_id = Column(String, primary_key=True)
    alternative_sequence = Column(String)
    chromosome = Column(String)
    data_sources = Column(String)
    xrefs = Column(ARRAY(String))
    position = Column(String)
    reference_sequence = Column(String)
    type = Column(String)
    variant_type = Column(String)
    observation_source = Column(String)


### Association tables ###

# Connections between cohort observations and external knowledge graph entities

class CohortReferencesPhenotype(Base):
    __tablename__ = 'cohort_references_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    cohort_id = Column(String, ForeignKey('cohort_phenotype.cohort_id'))
    hpo_id = Column(String, ForeignKey('phenotype.hpo_id'))


class CohortReferencesDisease(Base):
    __tablename__ = 'cohort_references_disease'
    id = Column(Integer, primary_key=True, autoincrement=True)
    cohort_id = Column(String, ForeignKey('cohort_phenotype.cohort_id'))
    mondo_id = Column(String, ForeignKey('disorder.mondo_id'))


class CohortReferencesProtein(Base):
    __tablename__ = 'cohort_references_protein'
    id = Column(Integer, primary_key=True, autoincrement=True)
    cohort_id = Column(String, ForeignKey('cohort_protein.cohort_id'))
    uniprot_id = Column(String, ForeignKey('protein.uniprot_id'))


class CohortReferencesMetabolite(Base):
    __tablename__ = 'cohort_references_metabolite'
    id = Column(Integer, primary_key=True, autoincrement=True)
    cohort_id = Column(String, ForeignKey('cohort_metabolite.cohort_id'))
    hmdb_id = Column(String, ForeignKey('metabolite.hmdb_id'))


class CohortReferencesVariant(Base):
    __tablename__ = 'cohort_references_variant'
    id = Column(Integer, primary_key=True, autoincrement=True)
    cohort_id = Column(String, ForeignKey('cohort_variant.cohort_id'))
    clinvar_id = Column(String, ForeignKey('genomic_variant.clinvar_id'))


# Calculated effects of cohort observations
class EffectsVariantProtein(Base):
    __tablename__ = 'edges_variant_protein'
    id = Column(Integer, primary_key=True, autoincrement=True)
    variant_id = Column(String, ForeignKey('cohort_variant.cohort_id'))
    protein_id = Column(String, ForeignKey('cohort_protein.cohort_id'))

    gwas_p_unadjusted = Column(Float)
    gwas_p_bonferroni = Column(Float)
    gwas_e_unspecified = Column(Float)


class EffectsVariantMetabolite(Base):
    __tablename__ = 'edges_variant_metabolite'
    id = Column(Integer, primary_key=True, autoincrement=True)
    variant_id = Column(String, ForeignKey('cohort_variant.cohort_id'))
    metabolite_id = Column(String, ForeignKey('cohort_metabolite.cohort_id'))

    gwas_p_unadjusted = Column(Float)
    gwas_p_bonferroni = Column(Float)
    gwas_e_unspecified = Column(Float)


class EffectsVariantPhenotype(Base):
    __tablename__ = 'edges_variant_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    variant_id = Column(String, ForeignKey('cohort_variant.cohort_id'))
    phenotype_id = Column(String, ForeignKey('cohort_phenotype.cohort_id'))

    gwas_p_unadjusted = Column(Float)
    gwas_p_bonferroni = Column(Float)
    gwas_e_odds_ratio = Column(Float)
    gwas_e_unspecified = Column(Float)


class EffectsProteinProtein(Base):
    __tablename__ = 'edges_protein_protein'
    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id_1 = Column(String, ForeignKey('cohort_protein.cohort_id'))
    protein_id_2 = Column(String, ForeignKey('cohort_protein.cohort_id'))

    pearson_p_unadjusted = Column(Float)
    pearson_p_bonferroni = Column(Float)
    pearson_p_benjamini_hb = Column(Float)
    pearson_p_benjamini_yek = Column(Float)
    pearson_e_r2 = Column(Float)

    spearman_p_unadjusted = Column(Float)
    spearman_p_bonferroni = Column(Float)
    spearman_p_benjamini_hb = Column(Float)
    spearman_p_benjamini_yek = Column(Float)
    spearman_e_rho = Column(Float)


class EffectsProteinMetabolite(Base):
    __tablename__ = 'edges_protein_metabolite'
    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(String, ForeignKey('cohort_protein.cohort_id'))
    metabolite_id = Column(String, ForeignKey('cohort_metabolite.cohort_id'))

    pearson_p_unadjusted = Column(Float)
    pearson_p_bonferroni = Column(Float)
    pearson_p_benjamini_hb = Column(Float)
    pearson_p_benjamini_yek = Column(Float)
    pearson_e_r2 = Column(Float)

    spearman_p_unadjusted = Column(Float)
    spearman_p_bonferroni = Column(Float)
    spearman_p_benjamini_hb = Column(Float)
    spearman_p_benjamini_yek = Column(Float)
    spearman_e_rho = Column(Float)


class EffectsProteinPhenotype(Base):
    __tablename__ = 'edges_protein_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    protein_id = Column(String, ForeignKey('cohort_protein.cohort_id'))
    phenotype_id = Column(String, ForeignKey('cohort_phenotype.cohort_id'))

    ttest_p_unadjusted = Column(Float)
    ttest_p_bonferroni = Column(Float)
    ttest_p_benjamini_hb = Column(Float)
    ttest_p_benjamini_yek = Column(Float)
    ttest_e_cohens_d = Column(Float)

    anova_p_unadjusted = Column(Float)
    anova_p_bonferroni = Column(Float)
    anova_p_benjamini_hb = Column(Float)
    anova_p_benjamini_yek = Column(Float)
    anova_e_np2 = Column(Float)

    mwu_p_unadjusted = Column(Float)
    mwu_p_bonferroni = Column(Float)
    mwu_p_benjamini_hb = Column(Float)
    mwu_p_benjamini_yek = Column(Float)
    mwu_e_r = Column(Float)

    kruskal_p_unadjusted = Column(Float)
    kruskal_p_bonferroni = Column(Float)
    kruskal_p_benjamini_hb = Column(Float)
    kruskal_p_benjamini_yek = Column(Float)
    kruskal_e_eta2 = Column(Float)

    pearson_p_unadjusted = Column(Float)
    pearson_p_bonferroni = Column(Float)
    pearson_p_benjamini_hb = Column(Float)
    pearson_p_benjamini_yek = Column(Float)
    pearson_e_r2 = Column(Float)

    spearman_p_unadjusted = Column(Float)
    spearman_p_bonferroni = Column(Float)
    spearman_p_benjamini_hb = Column(Float)
    spearman_p_benjamini_yek = Column(Float)
    spearman_e_rho = Column(Float)


class EffectsMetaboliteMetabolite(Base):
    __tablename__ = 'edges_metabolite_metabolite'
    id = Column(Integer, primary_key=True, autoincrement=True)
    metabolite_id_1 = Column(String, ForeignKey('cohort_metabolite.cohort_id'))
    metabolite_id_2 = Column(String, ForeignKey('cohort_metabolite.cohort_id'))

    pearson_p_unadjusted = Column(Float)
    pearson_p_bonferroni = Column(Float)
    pearson_p_benjamini_hb = Column(Float)
    pearson_p_benjamini_yek = Column(Float)
    pearson_e_r2 = Column(Float)

    spearman_p_unadjusted = Column(Float)
    spearman_p_bonferroni = Column(Float)
    spearman_p_benjamini_hb = Column(Float)
    spearman_p_benjamini_yek = Column(Float)
    spearman_e_rho = Column(Float)


class EffectsMetabolitePhenotype(Base):
    __tablename__ = 'edges_metabolite_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    metabolite_id = Column(String, ForeignKey('cohort_metabolite.cohort_id'))
    phenotype_id = Column(String, ForeignKey('cohort_phenotype.cohort_id'))

    ttest_p_unadjusted = Column(Float)
    ttest_p_bonferroni = Column(Float)
    ttest_p_benjamini_hb = Column(Float)
    ttest_p_benjamini_yek = Column(Float)
    ttest_e_cohens_d = Column(Float)

    anova_p_unadjusted = Column(Float)
    anova_p_bonferroni = Column(Float)
    anova_p_benjamini_hb = Column(Float)
    anova_p_benjamini_yek = Column(Float)
    anova_e_np2 = Column(Float)

    mwu_p_unadjusted = Column(Float)
    mwu_p_bonferroni = Column(Float)
    mwu_p_benjamini_hb = Column(Float)
    mwu_p_benjamini_yek = Column(Float)
    mwu_e_r = Column(Float)

    kruskal_p_unadjusted = Column(Float)
    kruskal_p_bonferroni = Column(Float)
    kruskal_p_benjamini_hb = Column(Float)
    kruskal_p_benjamini_yek = Column(Float)
    kruskal_e_eta2 = Column(Float)

    pearson_p_unadjusted = Column(Float)
    pearson_p_bonferroni = Column(Float)
    pearson_p_benjamini_hb = Column(Float)
    pearson_p_benjamini_yek = Column(Float)
    pearson_e_r2 = Column(Float)

    spearman_p_unadjusted = Column(Float)
    spearman_p_bonferroni = Column(Float)
    spearman_p_benjamini_hb = Column(Float)
    spearman_p_benjamini_yek = Column(Float)
    spearman_e_rho = Column(Float)


class EffectsPhenotypePhenotype(Base):
    __tablename__ = 'edges_phenotype_phenotype'
    id = Column(Integer, primary_key=True, autoincrement=True)
    phenotype_id_1 = Column(String, ForeignKey('cohort_phenotype.cohort_id'))
    phenotype_id_2 = Column(String, ForeignKey('cohort_phenotype.cohort_id'))

    chi2_p_unadjusted = Column(Float)
    chi2_p_bonferroni = Column(Float)
    chi2_p_benjamini_hb = Column(Float)
    chi2_p_benjamini_yek = Column(Float)
    chi2_e_cramers_v = Column(Float)
    chi2_e_phi = Column(Float)

    ttest_p_unadjusted = Column(Float)
    ttest_p_bonferroni = Column(Float)
    ttest_p_benjamini_hb = Column(Float)
    ttest_p_benjamini_yek = Column(Float)
    ttest_e_cohens_d = Column(Float)

    anova_p_unadjusted = Column(Float)
    anova_p_bonferroni = Column(Float)
    anova_p_benjamini_hb = Column(Float)
    anova_p_benjamini_yek = Column(Float)
    anova_e_np2 = Column(Float)

    mwu_p_unadjusted = Column(Float)
    mwu_p_bonferroni = Column(Float)
    mwu_p_benjamini_hb = Column(Float)
    mwu_p_benjamini_yek = Column(Float)
    mwu_e_r = Column(Float)

    kruskal_p_unadjusted = Column(Float)
    kruskal_p_bonferroni = Column(Float)
    kruskal_p_benjamini_hb = Column(Float)
    kruskal_p_benjamini_yek = Column(Float)
    kruskal_e_eta2 = Column(Float)

    pearson_p_unadjusted = Column(Float)
    pearson_p_bonferroni = Column(Float)
    pearson_p_benjamini_hb = Column(Float)
    pearson_p_benjamini_yek = Column(Float)
    pearson_e_r2 = Column(Float)

    spearman_p_unadjusted = Column(Float)
    spearman_p_bonferroni = Column(Float)
    spearman_p_benjamini_hb = Column(Float)
    spearman_p_benjamini_yek = Column(Float)
    spearman_e_rho = Column(Float)


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


class ProteinAssocGene(Base):
    __tablename__ = 'protein_associates_gene'
    id = Column(Integer, primary_key=True, autoincrement=True)
    uniprot_id = Column(String, ForeignKey('protein.uniprot_id'))
    entrez_id = Column(String, ForeignKey('gene.entrez_id'))


class VariantAssocGene(Base):
    __tablename__ = 'variant_associates_gene'
    id = Column(Integer, primary_key=True, autoincrement=True)
    clinvar_id = Column(String, ForeignKey('genomic_variant.clinvar_id'))
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


### New-style tables: single nodes table + parametric/nonparametric edge tables ###

class Node(Base):
    __tablename__ = 'nodes'
    node_id = Column(String, primary_key=True)
    display_name = Column(String)
    data_type = Column(String)
    node_group = Column(String, nullable=True)
    description = Column(String, nullable=True)
    xrefs = Column(String, nullable=True)


class EdgeParametric(Base):
    __tablename__ = 'edges_parametric'
    # node_id_1/node_id_2 are canonicalized (smaller node_id always first, see
    # insert_scores()) so this constraint catches a pair regardless of which
    # order it was loaded in -- an undirected edge is the same edge either way.
    __table_args__ = (
        UniqueConstraint('node_id_1', 'node_id_2', name='edges_parametric_pair_unique'),
    )
    id = Column(Integer, primary_key=True, autoincrement=True)
    node_id_1 = Column(String, ForeignKey('nodes.node_id'))
    node_id_2 = Column(String, ForeignKey('nodes.node_id'))
    p_value = Column(Float)
    effect_size = Column(Float, nullable=True)
    test_type = Column(String)


class EdgeNonparametric(Base):
    __tablename__ = 'edges_nonparametric'
    __table_args__ = (
        UniqueConstraint('node_id_1', 'node_id_2', name='edges_nonparametric_pair_unique'),
    )
    id = Column(Integer, primary_key=True, autoincrement=True)
    node_id_1 = Column(String, ForeignKey('nodes.node_id'))
    node_id_2 = Column(String, ForeignKey('nodes.node_id'))
    p_value = Column(Float)
    effect_size = Column(Float, nullable=True)
    test_type = Column(String)
