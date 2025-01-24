library(clusterProfiler)
library(enrichR)

table1 <- read.csv("C:\\Users\\piercetf\\Documents\\nparc_outputs_Oct2024.csv")
table2 <- read.csv("C:\\Users\\piercetf\\Documents\\nparc_unshared_fisetin_Oct2024.csv")

orep1 <- enrichGO(table2$PG.ProteinAccessions, 
         'org.Hs.eg.db',
         keyType='UNIPROT',
         ont='BP', 
         universe=table1$PG.ProteinAccessions)

dotplot(orep1)


orep1 <- enrichGO(table2$PG.ProteinAccessions, 
                  'org.Hs.eg.db',
                  keyType='UNIPROT',
                  ont='MF', 
                  universe=table1$PG.ProteinAccessions)

dotplot(orep1)

orep1 <- enrichGO(table2$PG.ProteinAccessions, 
                  'org.Hs.eg.db',
                  keyType='UNIPROT',
                  ont='CC', 
                  universe=table1$PG.ProteinAccessions)

dotplot(orep1)

dbs <- c("Pfam_InterPro_Domains", "Pfam_Domains_2019", "InterPro_Domains_2019",
         "Reactome_2022", "Reactome_Pathways_2024", "KEGG_2019_Human")

fis_enrich <- enrichr(table2$PG.Genes, dbs)

plotEnrich(fis_enrich[[1]], title="domains 1")

plotEnrich(fis_enrich[[2]], title="domains 2")

plotEnrich(fis_enrich[[3]], title="domains 3")

plotEnrich(fis_enrich[[4]], title="reactome 1")

plotEnrich(fis_enrich[[5]], title="reactome 2")

plotEnrich(fis_enrich[[6]], title="kegg")
