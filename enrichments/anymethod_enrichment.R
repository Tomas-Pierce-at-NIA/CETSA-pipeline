library(clusterProfiler)
library(enrichR)
library(org.Hs.eg.db)

background <- read.csv(r"(C:\Users\piercetf\Projects\senocetsa_outputs\nparc_Fisetin_Jan2025.csv)")
nparc <- read.csv(r"(C:\Users\piercetf\Projects\senocetsa_outputs\nparc_unshare_Fisetin_Jan2025.csv)")
ita <- read.csv(r"(C:\Users\piercetf\Projects\senocetsa_outputs\ITA_un_Fisetin_Jan2025.csv)")

any <- c(nparc$PG.Genes, ita$PG.Genes)

dbs <- c("Pfam_InterPro_Domains", "Pfam_Domains_2019", "InterPro_Domains_2019",
         "Reactome_2022", "Reactome_Pathways_2024", "KEGG_2019_Human")

fis_enrich <- enrichr(any, dbs, background=background$PG.Genes)


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



background <- read.csv(r"(C:\Users\piercetf\Projects\senocetsa_outputs\nparc_Quercetin_Jan2025.csv)")
nparc <- read.csv(r"(C:\Users\piercetf\Projects\senocetsa_outputs\nparc_unshare_Quercetin_Jan2025.csv)")
ita <- read.csv(r"(C:\Users\piercetf\Projects\senocetsa_outputs\ITA_un_Quercetin_Jan2025.csv)")

any <- c(nparc$PG.Genes, ita$PG.Genes)

dbs <- c("Pfam_InterPro_Domains", "Pfam_Domains_2019", "InterPro_Domains_2019",
         "Reactome_2022", "Reactome_Pathways_2024", "KEGG_2019_Human")

quer_enrich <- enrichr(any, dbs, background=background$PG.Genes)


plotEnrich(quer_enrich[[1]], 
           title="Pfam InterPro Quercetin",
           orderBy="Adjusted.P.value")

plotEnrich(quer_enrich[[2]], 
           title="Pfam Quercetin",
           orderBy="Adjusted.P.value")

plotEnrich(quer_enrich[[3]], 
           title="InterPro Quercetin",
           orderBy="Adjusted.P.value")

plotEnrich(quer_enrich[[4]], 
           title="reactome 1",
           orderBy="Adjusted.P.value")

plotEnrich(quer_enrich[[5]], 
           title="reactome 2",
           orderBy="Adjusted.P.value")

plotEnrich(quer_enrich[[6]], 
           title="kegg",
           orderBy="Adjusted.P.value")

table1_idents <- bitr(background$PG.ProteinAccessions, 
                      "UNIPROT",
                      "ENTREZID",
                      "org.Hs.eg.db")

anyprot <- c(nparc$PG.ProteinAccessions, ita$PG.ProteinAccessions)

table2_idents <- bitr(anyprot,
                      "UNIPROT",
                      "ENTREZID",
                      "org.Hs.eg.db")

# GO enrichment

orep1 <- enrichGO(table2_idents$ENTREZID, 
                  'org.Hs.eg.db',
                  ont='MF', 
                  universe=table1_idents$ENTREZID)

dotplot(orep1, title="Quercetin Molecular Function")

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

