SELECT *
FROM fisetin_nparc_detects
INNER JOIN Senescent_Diff_Expression AS S
ON "PG.Genes" = Gene
WHERE S.Pvalue < 0.05
ORDER BY "AVG Log2 Ratio"
