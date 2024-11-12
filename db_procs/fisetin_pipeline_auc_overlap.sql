SELECT A.Genes
FROM fisetin_auc_detects AS A
INNER JOIN fisetin_nparc_detects AS F
ON A.Genes = "PG.Genes"
UNION
SELECT A.Genes
FROM fisetin_auc_detects AS A
INNER JOIN fisetin_ita_detects AS F
ON A.Genes = "PG.Genes";