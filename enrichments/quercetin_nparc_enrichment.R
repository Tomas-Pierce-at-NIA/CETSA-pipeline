library(clusterProfiler)
library(enrichR)
library(org.Hs.eg.db)

table1 <- read.csv(r"(C:\Users\piercetf\Projects\example_output3\nparc_Quercetin_Jan2025.csv)")
table2 <- read.csv(r"(C:\Users\piercetf\Projects\example_output3\nparc_unshare_Quercetin_Jan2025.csv)")

table1_idents <- bitr(table1$PG.ProteinAccessions, 
                      "UNIPROT",
                      "ENTREZID",
                      "org.Hs.eg.db")

table2_idents <- bitr(table2$PG.ProteinAccessions,
                      "UNIPROT",
                      "ENTREZID",
                      "org.Hs.eg.db")

orep1 <- enrichGO(table2_idents$ENTREZID, 
                  'org.Hs.eg.db',
                  ont='MF', 
                  universe=table1_idents$ENTREZID)

dotplot(orep1, title="Quercetin MF")

orep2 <- enrichGO(table2_idents$ENTREZID, 
                  'org.Hs.eg.db',
                  ont='CC', 
                  universe=table1_idents$ENTREZID)

dotplot(orep2)

orep3 <- enrichGO(table2_idents$ENTREZID, 
                  'org.Hs.eg.db',
                  ont='BP', 
                  universe=table1_idents$ENTREZID)

dotplot(orep3)


dbs <- c("Pfam_InterPro_Domains", "Pfam_Domains_2019", "InterPro_Domains_2019",
         "Reactome_2022", "Reactome_Pathways_2024", "KEGG_2019_Human")

quer_enrich <- enrichr(table2$PG.Genes, dbs, background=table1$PG.Genes)

plotEnrich(quer_enrich[[1]], title="Pfam InterPro Quercetin")

plotEnrich(quer_enrich[[2]], title="Pfam Quercetin")

plotEnrich(quer_enrich[[3]], title="InterPro Quercetin")

plotEnrich(quer_enrich[[4]], title="reactome 1")

plotEnrich(quer_enrich[[5]], title="reactome 2")

plotEnrich(quer_enrich[[6]], title="kegg")

