SELECT "PG.Genes" AS GeneID, *
FROM fisetin_nparc_detects AS FN
FULL OUTER JOIN fisetin_ita_detects AS FI
USING ("PG.Genes")
FULL OUTER JOIN quercetin_nparc_detects AS QN
USING ("PG.Genes")
FULL OUTER JOIN quercetin_ita_detects AS QI
USING ("PG.Genes")
--FULL OUTER JOIN fisetin_auc_detects AS FA ON "PG.Genes" = FA.Genes
--FULL OUTER JOIN quercetin_auc_detects AS QA ON "PG.Genes" = QA.Genes
LEFT JOIN Senescent_Diff_Expression AS SDE
ON "PG.Genes" = SDE.Gene
ORDER BY GeneID
