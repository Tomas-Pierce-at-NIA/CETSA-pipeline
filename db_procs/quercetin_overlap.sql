CREATE VIEW Quercetin_Overlap AS SELECT "PG.Genes"
FROM quercetin_ita_detects
INTERSECT
SELECT "PG.Genes"
FROM quercetin_nparc_detects