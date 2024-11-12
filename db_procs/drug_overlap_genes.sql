CREATE VIEW Drug_Overlap AS WITH ita AS (SELECT "PG.Genes"
FROM fisetin_ita_detects
INTERSECT
SELECT "PG.Genes"
FROM quercetin_ita_detects),
nparc AS (
SELECT "PG.Genes"
FROM fisetin_nparc_detects
INTERSECT
SELECT "PG.Genes"
FROM quercetin_nparc_detects)
SELECT *
FROM nparc
FULL OUTER JOIN ita
USING ("PG.Genes")

