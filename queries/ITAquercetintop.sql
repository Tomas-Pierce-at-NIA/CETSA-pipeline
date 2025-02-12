SELECT *
FROM ITA_un_Quercetin_Jan2025 AS A
WHERE bh_pval = (
  SELECT min(bh_pval)
  FROM ITA_un_Quercetin_Jan2025 AS B
  WHERE A.PG_Genes = B.PG_Genes
  )
ORDER BY bh_pval
LIMIT 20