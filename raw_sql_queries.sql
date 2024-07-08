SELECT COUNT(genes.display_name)
FROM genes;


-- rename genome variant column to clinvar_id
ALTER TABLE genomic_variant
RENAME COLUMN "variant_primaryDomainId" TO clinvar_id;


-- this is for typeahead search
CREATE MATERIALIZED VIEW view_description_fts AS
SELECT 'disorder' AS source_table, mondo_id AS id, description, NULL AS display_name FROM disorders
UNION ALL
SELECT 'metabolite' AS source_table, hmdb_id AS id, description, display_name FROM metabolites
UNION ALL
SELECT 'gene' AS source_table, entrez_id AS id, description, display_name FROM genes
UNION ALL
SELECT 'protein' AS source_table, uniprot_id AS id, description, NULL AS display_name FROM proteins
UNION ALL
SELECT 'phenotype' AS source_table, hpo_id AS id, description, display_name FROM phenotypes;

-- drop the materialized view if it already exists
DROP MATERIALIZED VIEW view_description_fts;

SELECT *
FROM view_description_fts;

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


-- select all rows from the proteins table
SELECT *
FROM protein_associates_proteins
JOIN proteins ON protein_associates_proteins."uniprot_id_memberOne" = proteins.uniprot_id
WHERE "uniprot_id_memberOne" = 'uniprot.P31937' OR "uniprot_id_memberTwo" = 'uniprot.P31937';

-- check if a disorder has an OMIM xref
SELECT *
FROM disorders
WHERE EXISTS (
    SELECT 1
    FROM unnest(xrefs) AS xref
    WHERE xref LIKE 'omim.%'
) AND observation_source = 'CHRIS';

-- create an index on the uniprot_id column for faster lookups
CREATE INDEX idx_uniprot_id_1 ON effects_protein_protein(uniprot_id_1);
CREATE INDEX idx_uniprot_id_2 ON effects_protein_protein(uniprot_id_2);

-- select all rows from the effects_protein_protein table where a specific protein is involved
SELECT
    effects_protein_protein.*,
    proteins1.*,
    proteins2.*
FROM
    effects_protein_protein
JOIN
    proteins AS proteins1 ON effects_protein_protein."uniprot_id_1" = proteins1.uniprot_id
JOIN
    proteins AS proteins2 ON effects_protein_protein.uniprot_id_2 = proteins2.uniprot_id
WHERE
    effects_protein_protein.uniprot_id_1 = 'uniprot.P31937'
    OR effects_protein_protein.uniprot_id_2 = 'uniprot.P31937';

-- select all rows from the effects_phenotype_disorder table where a specific phenotype is involved
SELECT
    effects_phenotype_disorder.*,
    phenotypes.*,
    disorders.*
FROM effects_phenotype_disorder
JOIN phenotypes ON effects_phenotype_disorder.hpo_id = phenotypes.hpo_id
JOIN disorders ON effects_phenotype_disorder.mondo_id = disorders.mondo_id
WHERE effects_phenotype_disorder.hpo_id = 'hpo.0012398';


-- select all tables that contain a specific column. these  are the ones we have to go through to find the uniprot_id
SELECT table_name, column_name
FROM information_schema.columns
WHERE column_name LIKE 'uniprot_id%'
AND table_schema = 'public';
