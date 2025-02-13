library(clusterProfiler)
library(enrichR)
library(org.Hs.eg.db)

table2 <- read.csv(
  r"(C:\Users\piercetf\Projects\senocetsa_outputs\ITA_un_Fisetin_Jan2025.csv)"
  )

table1 <- read.csv(
  r"(C:\Users\piercetf\Projects\senocetsa_outputs\ITA_all_comps_Student_Oct2024.csv)"
)

dbs <- c("Pfam_InterPro_Domains", "Pfam_Domains_2019", "InterPro_Domains_2019",
         "Reactome_2022", "Reactome_Pathways_2024", "KEGG_2019_Human")

fis_enrich <- enrichr(table2$PG.Genes, dbs, background=unique(table1$PG.Genes))

plotEnrich(fis_enrich[[1]], 
           title="Pfam InterPro Fisetin",
           orderBy="Adjusted.P.value")

plotEnrich(fis_enrich[[2]], 
           title="Pfam Fisetin",
           orderBy="Adjusted.P.value")

plotEnrich(fis_enrich[[3]], 
           title="InterPro Fisetin",
           orderBy="Adjusted.P.value")

plotEnrich(fis_enrich[[4]], 
           title="reactome 1",
           orderBy="Adjusted.P.value")

plotEnrich(fis_enrich[[5]], 
           title="reactome 2",
           orderBy="Adjusted.P.value")

plotEnrich(fis_enrich[[6]], 
           title="kegg",
           orderBy="Adjusted.P.value")

# convert identities for GO enrichment

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

dotplot(orep1, title="Fisetin Molecular Function")

orep2 <- enrichGO(table2_idents$ENTREZID, 
                  'org.Hs.eg.db',
                  ont='CC', 
                  universe=table1_idents$ENTREZID)

dotplot(orep2, title="Fisetin Cellular Compartment")

orep3 <- enrichGO(table2_idents$ENTREZID, 
                  'org.Hs.eg.db',
                  ont='BP', 
                  universe=table1_idents$ENTREZID)

dotplot(orep3, title="Fisetin BP")



