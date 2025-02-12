SELECT *
FROM nparc_unshare_Fisetin_Jan2025 AS A
WHERE bh_pval = (
  SELECT min(bh_pval)
  FROM nparc_unshare_Fisetin_Jan2025 AS B
  WHERE A.PG_Genes = B.PG_Genes
  )
ORDER BY bh_pval, (T_infl_Treatment_1 - T_infl_Treatment_2)
LIMIT 20