CREATE VIEW Fisetin_any_detects AS SELECT DISTINCT "PG.Genes"
FROM fisetin_nparc_detects
FULL OUTER JOIN fisetin_ita_detects
USING ("PG.Genes");
CREATE VIEW Quercetin_any_detects AS SELECT DISTINCT "PG.Genes"
FROM quercetin_nparc_detects
FULL OUTER JOIN quercetin_ita_detects
USING ("PG.Genes");
