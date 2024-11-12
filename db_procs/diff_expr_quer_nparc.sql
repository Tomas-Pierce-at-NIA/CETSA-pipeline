SELECT *
FROM quercetin_nparc_detects
LEFT JOIN Senescent_Diff_Expression AS S
ON "PG.Genes" = Gene
WHERE S.Pvalue < 0.05
ORDER BY bh_pval