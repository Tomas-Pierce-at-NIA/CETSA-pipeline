CREATE VIEW Sen_Differential_Expression AS 
SELECT Gene, "AVG Log2 Ratio", Pvalue
FROM Senescent_Diff_Expression
WHERE "# Unique Total Peptides" > 1
AND Pvalue < 0.05