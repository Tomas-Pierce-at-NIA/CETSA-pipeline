SELECT A.Genes
FROM quercetin_auc_detects AS A
INNER JOIN quercetin_nparc_detects AS F
ON A.Genes = "PG.Genes"
UNION
SELECT A.Genes
FROM quercetin_auc_detects AS A
INNER JOIN quercetin_ita_detects AS F
ON A.Genes = "PG.Genes";