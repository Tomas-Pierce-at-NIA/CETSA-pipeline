SELECT FI."PG.Genes", 
       FI.Treatment, 
	   FI.Control, 
	   FI.bh_pval, 
	   --FI.log_fold_change_auc, 
	   LOG2(FI."fold change AUC vs control") AS log2_foldchange_auc,
	   QI."PG.Genes", 
	   QI.Treatment, 
	   QI.Control, 
	   QI.bh_pval, 
	   --QI.log_fold_change_auc,
	   LOG2(QI."fold change AUC vs control") AS log2_foldchange_auc,
	   FN."PG.Genes",
	   FN."Treatment 1",
	   FN."Treatment 2",
	   FN.bh_pval,
	   FN.T_infl_Treatment_1 - FN.T_infl_Treatment_2 AS scaled_delta_T_infl,
	   QN."PG.Genes",
	   QN."Treatment 1",
	   QN."Treatment 2",
	   QN.bh_pval,
	   QN.T_infl_Treatment_1 - QN.T_infl_Treatment_2 AS scaled_delta_T_infl,
	   SDE.Gene,
	   SDE.Qvalue,
	   SDE."# Unique Total Peptides",
	   SDE."AVG Log2 Ratio",
	   GO."FROM",
	   GO."Protein names",
	   GO."Gene Ontology (biological process)",
	   GO."Gene Ontology (cellular component)",
	   GO."Gene Ontology (molecular function)"
FROM fisetin_ita_detects AS FI
FULL OUTER JOIN quercetin_ita_detects AS QI
ON FI."PG.Genes" = QI."PG.Genes" AND FI.Control = QI.Control
FULL OUTER JOIN fisetin_nparc_detects AS FN
ON (FN."PG.Genes" = FI."PG.Genes" OR FN."PG.Genes" = QI."PG.Genes") 
    AND (FN."Treatment 2" = FI.Control OR FN."Treatment 2" = QI.Control)
FULL OUTER JOIN quercetin_nparc_detects AS QN
ON (QN."PG.Genes" = FI."PG.Genes" OR QN."PG.Genes" = QI."PG.Genes" OR QN."PG.Genes" = FN."PG.Genes") 
    AND (QN."Treatment 2" = FI.Control OR QN."Treatment 2" = QI.Control OR QN."Treatment 2" = FN."Treatment 2")
LEFT JOIN Senescent_Diff_Expression AS SDE
ON (SDE.Gene = FI."PG.Genes" OR SDE.GENE = QI."PG.Genes" OR SDE.Gene = FN."PG.Genes" OR SDE.Gene = QN."PG.Genes")
LEFT JOIN pipeline_hits_GO AS GO
ON (GO."FROM" = FI."PG.Genes" OR GO."FROM" = QI."PG.Genes" OR GO."FROM" = FN."PG.Genes" OR GO."FROM" = QN."PG.Genes")
ORDER BY FI."PG.Genes" IS NULL, 
         FI."PG.Genes", 
		 QI."PG.Genes" IS NULL,
		 QI."PG.Genes",
		 FN."PG.Genes" IS NULL,
		 FN."PG.Genes", 
		 QN."PG.Genes" IS NULL,
		 QN."PG.Genes"