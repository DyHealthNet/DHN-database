SELECT COUNT(gene.display_name)
FROM gene;


-- typeahead search using the cohort information
CREATE MATERIALIZED VIEW view_cohort_fts AS
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

-- create materialized view for references edges
CREATE VIEW view_references_edges AS
SELECT 'protein' AS source_table, cohort_id, uniprot_id as reference_id FROM cohort_references_protein
UNION ALL
SELECT 'metabolite' AS source_table, cohort_id, hmdb_id as reference_id FROM cohort_references_metabolite
UNION ALL
SELECT 'phenotype' AS source_table, cohort_id, hpo_id as reference_id FROM cohort_references_phenotype
UNION ALL
SELECT 'disease' AS source_table, cohort_id, mondo_id as reference_id FROM cohort_references_disease;

WITH nodes AS (
    SELECT UNNEST(ARRAY['uniprot.Q9BUT1', 'uniprot.Q3SXY7', 'hmdb.HMDB0000011', 'uniprot.P22087', 'uniprot.P28908',
        'uniprot.Q13421', 'uniprot.Q9UM07', 'uniprot.Q96DN0', 'hmdb.HMDB0008189']) AS node_id
)
SELECT * FROM view_references_edges
WHERE reference_id IN (SELECT node_id FROM nodes);

-- create view for association edges
CREATE VIEW view_associations_edges AS
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
SELECT genomic_variant AS source_id, entrez_id AS target_id FROM variant_affects_gene;


SELECT *
FROM view_associations_edges
WHERE source_id = 'uniprot.P31937' OR target_id = 'uniprot.P31937';



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

SET enable_indexscan = off;
SET enable_bitmapscan = off;

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

-- test another
SELECT "view_associations_edges"."source_id", "view_associations_edges"."target_id" FROM "view_associations_edges" WHERE ("view_associations_edges"."source_id" IN ('hmdb.HMDB0013129', 'hmdb.HMDB0007920', 'hmdb.HMDB0008013', 'hmdb.HMDB0008109', 'hmdb.HMDB0008234', 'hmdb.HMDB0007926', 'hmdb.HMDB0008173', 'hmdb.HMDB0008300', 'hmdb.HMDB0008076', 'hmdb.HMDB0008590', 'uniprot.Q4V9L6', 'hmdb.HMDB0008788', 'hmdb.HMDB0008170', 'hmdb.HMDB0008042', 'hmdb.HMDB0000756', 'hmdb.HMDB0013134', 'uniprot.P07327', 'hmdb.HMDB0008107', 'hmdb.HMDB0008039', 'hmdb.HMDB0006317', 'hmdb.HMDB0008135', 'hmdb.HMDB0007951', 'hmdb.HMDB0000904', 'hmdb.HMDB0008494', 'hmdb.HMDB0008014', 'hmdb.HMDB0008071', 'hmdb.HMDB0008462', 'hmdb.HMDB0002014', 'hmdb.HMDB0007924', 'hmdb.HMDB0008365', 'hmdb.HMDB0008138', 'hmdb.HMDB0008017', 'hmdb.HMDB0008073', 'hmdb.HMDB0007889', 'hmdb.HMDB0008690', 'hmdb.HMDB0008203', 'uniprot.Q16206', 'uniprot.P18031', 'hmdb.HMDB0008398', 'hmdb.HMDB0008142', 'hmdb.HMDB0008529', 'hmdb.HMDB0008044', 'hmdb.HMDB0008106', 'hmdb.HMDB0008299', 'hmdb.HMDB0007923', 'hpo.0000651', 'hmdb.HMDB0013329', 'hmdb.HMDB0008070', 'hmdb.HMDB0008496', 'hmdb.HMDB0008623', 'hmdb.HMDB0007894', 'uniprot.Q9BUJ0', 'hmdb.HMDB0000593', 'hmdb.HMDB0007982', 'hmdb.HMDB0008205', 'hmdb.HMDB0008561', 'uniprot.P15144', 'uniprot.Q13449', 'hmdb.HMDB0008268', 'hmdb.HMDB0007983', 'hmdb.HMDB0008011', 'hmdb.HMDB0008657', 'hmdb.HMDB0008722', 'hmdb.HMDB0013415', 'hmdb.HMDB0013207', 'hmdb.HMDB0008331', 'hmdb.HMDB0008237', 'hmdb.HMDB0007979', 'uniprot.Q9H4F8', 'uniprot.P29372', 'hmdb.HMDB0240588', 'mondo.0004247', 'hmdb.HMDB0008074', 'hmdb.HMDB0008172', 'hmdb.HMDB0008269', 'uniprot.P05783', 'hmdb.HMDB0008756', 'hmdb.HMDB0008559', 'hmdb.HMDB0007888', 'hmdb.HMDB0008202', 'hmdb.HMDB0008206', 'hmdb.HMDB0007986', 'hmdb.HMDB0007892', 'uniprot.Q9Y2Q3', 'hmdb.HMDB0008169', 'hmdb.HMDB0008103', 'hmdb.HMDB0008018', 'hmdb.HMDB0008429') AND "view_associations_edges"."target_id" IN ('hmdb.HMDB0013129', 'hmdb.HMDB0007920', 'hmdb.HMDB0008013', 'hmdb.HMDB0008109', 'hmdb.HMDB0008234', 'hmdb.HMDB0007926', 'hmdb.HMDB0008173', 'hmdb.HMDB0008300', 'hmdb.HMDB0008076', 'hmdb.HMDB0008590', 'uniprot.Q4V9L6', 'hmdb.HMDB0008788', 'hmdb.HMDB0008170', 'hmdb.HMDB0008042', 'hmdb.HMDB0000756', 'hmdb.HMDB0013134', 'uniprot.P07327', 'hmdb.HMDB0008107', 'hmdb.HMDB0008039', 'hmdb.HMDB0006317', 'hmdb.HMDB0008135', 'hmdb.HMDB0007951', 'hmdb.HMDB0000904', 'hmdb.HMDB0008494', 'hmdb.HMDB0008014', 'hmdb.HMDB0008071', 'hmdb.HMDB0008462', 'hmdb.HMDB0002014', 'hmdb.HMDB0007924', 'hmdb.HMDB0008365', 'hmdb.HMDB0008138', 'hmdb.HMDB0008017', 'hmdb.HMDB0008073', 'hmdb.HMDB0007889', 'hmdb.HMDB0008690', 'hmdb.HMDB0008203', 'uniprot.Q16206', 'uniprot.P18031', 'hmdb.HMDB0008398', 'hmdb.HMDB0008142', 'hmdb.HMDB0008529', 'hmdb.HMDB0008044', 'hmdb.HMDB0008106', 'hmdb.HMDB0008299', 'hmdb.HMDB0007923', 'hpo.0000651', 'hmdb.HMDB0013329', 'hmdb.HMDB0008070', 'hmdb.HMDB0008496', 'hmdb.HMDB0008623', 'hmdb.HMDB0007894', 'uniprot.Q9BUJ0', 'hmdb.HMDB0000593', 'hmdb.HMDB0007982', 'hmdb.HMDB0008205', 'hmdb.HMDB0008561', 'uniprot.P15144', 'uniprot.Q13449', 'hmdb.HMDB0008268', 'hmdb.HMDB0007983', 'hmdb.HMDB0008011', 'hmdb.HMDB0008657', 'hmdb.HMDB0008722', 'hmdb.HMDB0013415', 'hmdb.HMDB0013207', 'hmdb.HMDB0008331', 'hmdb.HMDB0008237', 'hmdb.HMDB0007979', 'uniprot.Q9H4F8', 'uniprot.P29372', 'hmdb.HMDB0240588', 'mondo.0004247', 'hmdb.HMDB0008074', 'hmdb.HMDB0008172', 'hmdb.HMDB0008269', 'uniprot.P05783', 'hmdb.HMDB0008756', 'hmdb.HMDB0008559', 'hmdb.HMDB0007888', 'hmdb.HMDB0008202', 'hmdb.HMDB0008206', 'hmdb.HMDB0007986', 'hmdb.HMDB0007892', 'uniprot.Q9Y2Q3', 'hmdb.HMDB0008169', 'hmdb.HMDB0008103', 'hmdb.HMDB0008018', 'hmdb.HMDB0008429')); args=('hmdb.HMDB0013129', 'hmdb.HMDB0007920', 'hmdb.HMDB0008013', 'hmdb.HMDB0008109', 'hmdb.HMDB0008234', 'hmdb.HMDB0007926', 'hmdb.HMDB0008173', 'hmdb.HMDB0008300', 'hmdb.HMDB0008076', 'hmdb.HMDB0008590', 'uniprot.Q4V9L6', 'hmdb.HMDB0008788', 'hmdb.HMDB0008170', 'hmdb.HMDB0008042', 'hmdb.HMDB0000756', 'hmdb.HMDB0013134', 'uniprot.P07327', 'hmdb.HMDB0008107', 'hmdb.HMDB0008039', 'hmdb.HMDB0006317', 'hmdb.HMDB0008135', 'hmdb.HMDB0007951', 'hmdb.HMDB0000904', 'hmdb.HMDB0008494', 'hmdb.HMDB0008014', 'hmdb.HMDB0008071', 'hmdb.HMDB0008462', 'hmdb.HMDB0002014', 'hmdb.HMDB0007924', 'hmdb.HMDB0008365', 'hmdb.HMDB0008138', 'hmdb.HMDB0008017', 'hmdb.HMDB0008073', 'hmdb.HMDB0007889', 'hmdb.HMDB0008690', 'hmdb.HMDB0008203', 'uniprot.Q16206', 'uniprot.P18031', 'hmdb.HMDB0008398', 'hmdb.HMDB0008142', 'hmdb.HMDB0008529', 'hmdb.HMDB0008044', 'hmdb.HMDB0008106', 'hmdb.HMDB0008299', 'hmdb.HMDB0007923', 'hpo.0000651', 'hmdb.HMDB0013329', 'hmdb.HMDB0008070', 'hmdb.HMDB0008496', 'hmdb.HMDB0008623', 'hmdb.HMDB0007894', 'uniprot.Q9BUJ0', 'hmdb.HMDB0000593', 'hmdb.HMDB0007982', 'hmdb.HMDB0008205', 'hmdb.HMDB0008561', 'uniprot.P15144', 'uniprot.Q13449', 'hmdb.HMDB0008268', 'hmdb.HMDB0007983', 'hmdb.HMDB0008011', 'hmdb.HMDB0008657', 'hmdb.HMDB0008722', 'hmdb.HMDB0013415', 'hmdb.HMDB0013207', 'hmdb.HMDB0008331', 'hmdb.HMDB0008237', 'hmdb.HMDB0007979', 'uniprot.Q9H4F8', 'uniprot.P29372', 'hmdb.HMDB0240588', 'mondo.0004247', 'hmdb.HMDB0008074', 'hmdb.HMDB0008172', 'hmdb.HMDB0008269', 'uniprot.P05783', 'hmdb.HMDB0008756', 'hmdb.HMDB0008559', 'hmdb.HMDB0007888', 'hmdb.HMDB0008202', 'hmdb.HMDB0008206', 'hmdb.HMDB0007986', 'hmdb.HMDB0007892', 'uniprot.Q9Y2Q3', 'hmdb.HMDB0008169', 'hmdb.HMDB0008103', 'hmdb.HMDB0008018', 'hmdb.HMDB0008429', 'hmdb.HMDB0013129', 'hmdb.HMDB0007920', 'hmdb.HMDB0008013', 'hmdb.HMDB0008109', 'hmdb.HMDB0008234', 'hmdb.HMDB0007926', 'hmdb.HMDB0008173', 'hmdb.HMDB0008300', 'hmdb.HMDB0008076', 'hmdb.HMDB0008590', 'uniprot.Q4V9L6', 'hmdb.HMDB0008788', 'hmdb.HMDB0008170', 'hmdb.HMDB0008042', 'hmdb.HMDB0000756', 'hmdb.HMDB0013134', 'uniprot.P07327', 'hmdb.HMDB0008107', 'hmdb.HMDB0008039', 'hmdb.HMDB0006317', 'hmdb.HMDB0008135', 'hmdb.HMDB0007951', 'hmdb.HMDB0000904', 'hmdb.HMDB0008494', 'hmdb.HMDB0008014', 'hmdb.HMDB0008071', 'hmdb.HMDB0008462', 'hmdb.HMDB0002014', 'hmdb.HMDB0007924', 'hmdb.HMDB0008365', 'hmdb.HMDB0008138', 'hmdb.HMDB0008017', 'hmdb.HMDB0008073', 'hmdb.HMDB0007889', 'hmdb.HMDB0008690', 'hmdb.HMDB0008203', 'uniprot.Q16206', 'uniprot.P18031', 'hmdb.HMDB0008398', 'hmdb.HMDB0008142', 'hmdb.HMDB0008529', 'hmdb.HMDB0008044', 'hmdb.HMDB0008106', 'hmdb.HMDB0008299', 'hmdb.HMDB0007923', 'hpo.0000651', 'hmdb.HMDB0013329', 'hmdb.HMDB0008070', 'hmdb.HMDB0008496', 'hmdb.HMDB0008623', 'hmdb.HMDB0007894', 'uniprot.Q9BUJ0', 'hmdb.HMDB0000593', 'hmdb.HMDB0007982', 'hmdb.HMDB0008205', 'hmdb.HMDB0008561', 'uniprot.P15144', 'uniprot.Q13449', 'hmdb.HMDB0008268', 'hmdb.HMDB0007983', 'hmdb.HMDB0008011', 'hmdb.HMDB0008657', 'hmdb.HMDB0008722', 'hmdb.HMDB0013415', 'hmdb.HMDB0013207', 'hmdb.HMDB0008331', 'hmdb.HMDB0008237', 'hmdb.HMDB0007979', 'uniprot.Q9H4F8', 'uniprot.P29372', 'hmdb.HMDB0240588', 'mondo.0004247', 'hmdb.HMDB0008074', 'hmdb.HMDB0008172', 'hmdb.HMDB0008269', 'uniprot.P05783', 'hmdb.HMDB0008756', 'hmdb.HMDB0008559', 'hmdb.HMDB0007888', 'hmdb.HMDB0008202', 'hmdb.HMDB0008206', 'hmdb.HMDB0007986', 'hmdb.HMDB0007892', 'uniprot.Q9Y2Q3', 'hmdb.HMDB0008169', 'hmdb.HMDB0008103', 'hmdb.HMDB0008018', 'hmdb.HMDB0008429');