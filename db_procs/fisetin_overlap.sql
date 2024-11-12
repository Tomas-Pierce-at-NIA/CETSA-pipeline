CREATE VIEW Fisetin_Overlaps AS SELECT "PG.Genes"
FROM fisetin_ita_detects
INTERSECT
SELECT "PG.Genes"
FROM fisetin_nparc_detects
