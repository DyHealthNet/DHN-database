from sqlalchemy import text, Table
from utils.models import *


def add_views(session):
    # check if the view already exists and if so, update it
    view_exists = session.execute(text("SELECT to_regclass('view_description_fts')")).scalar()
    if view_exists is not None:
        session.execute(text("REFRESH MATERIALIZED VIEW view_description_fts;"))
        session.commit()
        print("View view_description_fts already exists. Refreshed.")
    else:
        # sql alchemy doesn't support creating views, so we have to use raw sql
        view_sql = """
        CREATE MATERIALIZED VIEW view_description_fts AS
        SELECT 'cohort_protein' AS source_table, cohort_id AS id, description, 
                display_name, xrefs FROM cohort_protein
        UNION ALL
        SELECT 'cohort_metabolite' AS source_table, cohort_id AS id, description, 
                display_name, xrefs FROM cohort_metabolite
        UNION ALL
        SELECT 'cohort_phenotype' AS source_table, cohort_id AS id, description, 
                display_name, xrefs FROM cohort_phenotype
        UNION ALL
        SELECT 'cohort_variant' AS source_table, cohort_id AS id, description, 
                display_name, xrefs FROM cohort_variant;
        """
        session.execute(text(view_sql))
        print("Created view view_description_fts.")

    # create new view called view_references_edges
    view_exists = session.execute(text("SELECT to_regclass('view_references_edges')")).scalar()
    if view_exists is None:
        view_sql = """
        CREATE VIEW view_references_edges AS
        SELECT 'protein' AS source_table, cohort_id, uniprot_id as reference_id FROM cohort_references_protein
        UNION ALL
        SELECT 'metabolite' AS source_table, cohort_id, hmdb_id as reference_id FROM cohort_references_metabolite
        UNION ALL
        SELECT 'phenotype' AS source_table, cohort_id, hpo_id as reference_id FROM cohort_references_phenotype
        UNION ALL
        SELECT 'disease' AS source_table, cohort_id, mondo_id as reference_id FROM cohort_references_disease
        UNION ALL
        SELECT 'variant' AS source_table, cohort_id, clinvar_id as reference_id FROM cohort_references_variant;
        """
        session.execute(text(view_sql))
        print("Created view view_references_edges.")

    # create new view called external_node_ids
    view_exists = session.execute(text("SELECT to_regclass('view_external_nodes')")).scalar()
    if view_exists is None:
        view_sql = """
            CREATE VIEW view_external_nodes AS
            SELECT mondo_id as "node_id", 'disorder' as source_table FROM disorder
            UNION ALL
            SELECT entrez_id as "node_id", 'gene' as source_table FROM gene
            UNION ALL
            SELECT hmdb_id as "node_id", 'metabolite' as source_table FROM metabolite
            UNION ALL
            SELECT hpo_id as "node_id", 'phenotype' as source_table FROM phenotype
            UNION ALL
            SELECT uniprot_id as "node_id", 'protein' as source_table FROM protein
            UNION ALL
            SELECT clinvar_id as "node_id", 'genomic_variant' as source_table FROM genomic_variant;
        """
        session.execute(text(view_sql))
        print("Created view external_node_ids")

    # create new view called view_associations_edges
    view_exists = session.execute(text("SELECT to_regclass('view_associations_edges')")).scalar()
    if view_exists is None:
        view_sql = """
        CREATE MATERIALIZED VIEW view_associations_edges AS
        SELECT uniprot_id_1 AS source_id, uniprot_id_2 AS target_id FROM protein_associates_protein
        UNION ALL
        SELECT uniprot_id AS source_id, hmdb_id AS target_id FROM protein_associates_metabolite
        UNION ALL
        SELECT mondo_id AS source_id, hpo_id AS target_id FROM disorder_associates_phenotype
        UNION ALL
        SELECT entrez_id AS source_id, mondo_id AS target_id FROM gene_associates_disorder
        UNION ALL
        SELECT hmdb_id AS source_id, mondo_id AS target_id FROM metabolite_associates_disorder
        UNION ALL
        SELECT clinvar_id AS source_id, entrez_id AS target_id FROM variant_associates_gene;
        """
        session.execute(text(view_sql))
        print("Created view view_associations_edges.")
    else:
        session.execute(text("REFRESH MATERIALIZED VIEW view_description_fts;"))
        print("View view_associations_edges already exists. Refreshed.")
    session.commit()


def add_indexes(session, engine, metadata):
    # Protein indexes for quick search
    idx_uniprot_id_1 = Index('idx_uniprot_id_1', EffectsProteinProtein.protein_id_1)
    # check if the index already exists
    if not session.execute(text("SELECT to_regclass('idx_uniprot_id_1')")).scalar():
        idx_uniprot_id_1.create(engine)

    idx_uniprot_id_2 = Index('idx_uniprot_id_2', EffectsProteinProtein.protein_id_2)
    if not session.execute(text("SELECT to_regclass('idx_uniprot_id_2')")).scalar():
        idx_uniprot_id_2.create(engine)

    idx_effects_protein_pheno = Index('idx_effects_protein_pheno', EffectsProteinPhenotype.protein_id)
    if not session.execute(text("SELECT to_regclass('idx_effects_protein_pheno')")).scalar():
        idx_effects_protein_pheno.create(engine)

    idx_effects_protein_metabo = Index('idx_effects_protein_metabo', EffectsProteinMetabolite.protein_id)
    if not session.execute(text("SELECT to_regclass('idx_effects_protein_metabo')")).scalar():
        idx_effects_protein_metabo.create(engine)

    # Index for quick typeahead search
    view_description_fts = Table('view_description_fts', metadata, autoload_with=engine)
    idx_display_name_fts = Index('idx_display_name_fts', view_description_fts.c.display_name)
    if not session.execute(text("SELECT to_regclass('idx_display_name_fts')")).scalar():
        idx_display_name_fts.create(engine)

    # add index for view_associations_edges
    view_associations_edges = Table('view_associations_edges', metadata, autoload_with=engine)
    idx_assoc_source_id = Index('idx_source_id', view_associations_edges.c.source_id)
    if not session.execute(text("SELECT to_regclass('idx_source_id')")).scalar():
        idx_assoc_source_id.create(engine)

    idx_assoc_target_id = Index('idx_target_id', view_associations_edges.c.target_id)
    if not session.execute(text("SELECT to_regclass('idx_target_id')")).scalar():
        idx_assoc_target_id.create(engine)

    # add the last index that doesn't work well with sqlalchemy
    if session.execute(text("SELECT to_regclass('idx_description_fts')")).scalar():
        print("Index idx_description_fts already exists.")
        return
    session.execute(text("CREATE INDEX idx_description_fts "
                         "ON view_description_fts USING gin(to_tsvector('english', description));"))
    print("Created indexes")
    session.commit()
