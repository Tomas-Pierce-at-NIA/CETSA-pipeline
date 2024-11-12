WITH Fisetin_Hits AS (SELECT "PG.Genes"
	FROM fisetin_nparc_detects UNION
	SELECT "PG.Genes"
	FROM fisetin_ita_detects
	),
	Quercetin_Hits AS (SELECT "PG.Genes"
	FROM quercetin_nparc_detects UNION
	SELECT "PG.Genes"
	FROM quercetin_ita_detects
	)
SELECT "PG.Genes"
FROM Fisetin_Hits
INTERSECT SELECT "PG.Genes"
FROM Quercetin_Hits;

WITH Fisetin_Hits AS (SELECT "PG.Genes"
	FROM fisetin_nparc_detects
	),
	Quercetin_Hits AS (SELECT "PG.Genes"
	FROM quercetin_nparc_detects
	)
SELECT "PG.Genes"
FROM Fisetin_Hits
INTERSECT SELECT "PG.Genes"
FROM Quercetin_Hits;

WITH Fisetin_Hits AS (SELECT "PG.Genes"
	FROM fisetin_ita_detects
	),
	Quercetin_Hits AS (SELECT "PG.Genes"
	FROM quercetin_ita_detects
	)
SELECT "PG.Genes"
FROM Fisetin_Hits
INTERSECT SELECT "PG.Genes"
FROM Quercetin_Hits;
