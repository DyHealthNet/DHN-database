SELECT COUNT(gene.display_name)
FROM gene;


-- rename genome variant column to clinvar_id
ALTER TABLE genomic_variant
RENAME COLUMN "variant_primaryDomainId" TO clinvar_id;


-- typeahead search using the cohort information
CREATE VIEW view_cohort_fts AS
SELECT 'cohort_protein' AS source_table, cohort_id AS id, description, display_name FROM cohort_protein
UNION ALL
SELECT 'cohort_metabolite' AS source_table, cohort_id AS id, description, display_name FROM cohort_metabolite
UNION ALL
SELECT 'cohort_phenotype' AS source_table, cohort_id AS id, description, display_name FROM cohort_phenotype;

SELECT *
FROM view_cohort_fts
WHERE view_cohort_fts.display_name IS NULL;

DROP VIEW view_cohort_fts;

CREATE INDEX idx_description_fts ON view_description_fts USING gin(to_tsvector('english', description));
CREATE INDEX idx_display_name ON view_description_fts (display_name);

-- list all indices
SELECT *
FROM pg_indexes
WHERE schemaname = 'public';


-- also for typeahead search
SELECT *
FROM view_description_fts
WHERE to_tsvector('english', description) @@ plainto_tsquery('english', 'low body');

SELECT *
FROM view_description_fts
WHERE display_name ILIKE 'ornithine%';


-- select all rows from the proteins table
SELECT *
FROM protein_associates_protein
JOIN protein ON protein_associates_protein.uniprot_id_1 = protein.uniprot_id
WHERE uniprot_id_1 = 'uniprot.P31937' OR uniprot_id_2 = 'uniprot.P31937';

-- check if a disorder has an OMIM xref
SELECT *
FROM disorder
WHERE EXISTS (
    SELECT 1
    FROM unnest(xrefs) AS xref
    WHERE xref LIKE 'omim.%'
) AND observation_source = 'CHRIS';

-- create an index on the uniprot_id column for faster lookups
CREATE INDEX idx_uniprot_id_1 ON effects_protein_protein(protein_id_1);
CREATE INDEX idx_uniprot_id_2 ON effects_protein_protein(protein_id_2);

-- select all rows from the effects_protein_protein table where a specific protein is involved
EXPLAIN SELECT
    effects_protein_protein.*,
    proteins1.*,
    proteins2.*
FROM
    effects_protein_protein
JOIN
    protein AS proteins1 ON effects_protein_protein.protein_id_1 = proteins1.uniprot_id
JOIN
    protein AS proteins2 ON effects_protein_protein.protein_id_2 = proteins2.uniprot_id
WHERE
    effects_protein_protein.protein_id_1 = 'uniprot.P31937'
    OR effects_protein_protein.protein_id_2 = 'uniprot.P31937';

-- select all rows from the effects_phenotype_disorder table where a specific phenotype is involved
SELECT
    effects_phenotype_phenotype.*,
    phenotype.*,
    disorder.*
FROM effects_phenotype_phenotype
JOIN phenotype ON effects_phenotype_phenotype.phenotype_id_1 = phenotype.hpo_id
JOIN disorder ON effects_phenotype_phenotype.phenotype_id_2 = disorder.mondo_id
WHERE effects_phenotype_phenotype.phenotype_id_1 = 'hpo.0012398';


-- select all tables that contain a specific column. these  are the ones we have to go through to find the uniprot_id
SELECT table_name, column_name
FROM information_schema.columns
WHERE column_name LIKE 'uniprot_id%'
AND table_schema = 'public';


WITH nodes AS (
    SELECT UNNEST(ARRAY['uniprot.Q9BUT1', 'uniprot.Q3SXY7', 'hmdb.HMDB0000011', 'uniprot.P22087', 'uniprot.P28908',
        'uniprot.Q13421', 'uniprot.Q9UM07', 'uniprot.Q96DN0', 'hmdb.HMDB0008189']) AS node_id
)
SELECT e0.uniprot_id AS source,
       e0.hmdb_id AS target
FROM protein_associates_metabolite e0
JOIN nodes n1 ON e0.uniprot_id = n1.node_id
JOIN nodes n2 ON e0.hmdb_id = n2.node_id

UNION ALL

SELECT e1.uniprot_id_1 AS source,
       e1.uniprot_id_2 AS target
FROM protein_associates_protein e1
JOIN nodes n1 ON e1.uniprot_id_1 = n1.node_id
JOIN nodes n2 ON e1.uniprot_id_2 = n2.node_id;

SELECT *
FROM view_description_fts
WHERE to_tsvector('english', description) @@ plainto_tsquery('english','A2 receptor') OR
id ILIKE 'A2 receptor%' OR
display_name ILIKE 'A2 receptor%'
LIMIT 10;

-- some more indexes
CREATE INDEX idx_effects_protein_pheno ON effects_protein_phenotype(protein_id);
CREATE INDEX idx_effects_protein_metab ON effects_protein_metabolite(protein_id);


-- test queries

SELECT *
FROM "effects_protein_protein"
WHERE ("effects_protein_protein"."protein_id_1" = 'x0so1193' OR
       "effects_protein_protein"."protein_id_2" = 'x0so1193')
ORDER BY "effects_protein_protein"."p_value" ASC
LIMIT 10;

SELECT *
FROM "effects_protein_metabolite"
WHERE "effects_protein_metabolite"."protein_id" = 'x0so1193'
ORDER BY "effects_protein_metabolite"."p_value" ASC LIMIT 10;

SELECT *
FROM "effects_protein_phenotype"
WHERE "effects_protein_phenotype"."protein_id" = 'x0so1193'
ORDER BY "effects_protein_phenotype"."p_value" ASC
LIMIT 10;