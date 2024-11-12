
library(enrichR)

fisetin_name <- "C:/Users/piercetf/Documents/nparc_unshared_fisetin_Oct2024.csv"
quercetin_name <- "C:/Users/piercetf/Documents/nparc_unshared_quercetin_Oct2024.csv"
fisetin_name2 <- "C:/Users/piercetf/Documents/ITA_fisetin_Student_Oct2024.csv"
quercetin_name2 <- "C:/Users/piercetf/Documents/ITA_quercetin_Student_Oct2024.csv"

nparc_fisetin_data <- read.csv(fisetin_name)
nparc_quercetin_data <- read.csv(quercetin_name)
ita_fisetin_data <- read.csv(fisetin_name2)
ita_quercetin_data <- read.csv(quercetin_name2)

dbs <- listEnrichrDbs()

checkdbs <- c("GO_Biological_Process_2023",
              "GO_Cellular_Component_2023",
              "GO_Molecular_Function_2023",
              "Reactome_2022",
              "KEGG_2019_Human",
              "InterPro_Domains_2019")

nf <- enrichr(nparc_fisetin_data$PG.Genes, checkdbs)

pdf("C:/Users/piercetf/Documents/enrichr_fisetin_nparc.pdf")

plotEnrich(nf[[1]], title="fisetin BP N")
plotEnrich(nf[[2]], title="fisetin CC N")
plotEnrich(nf[[3]], title="fisetin MF N")
plotEnrich(nf[[4]], title="fisetin R-ome N")
plotEnrich(nf[[5]], title="fisetin kegg N")
plotEnrich(nf[[6]], title="fisetin InterPro N")

dev.off()

fi <- enrichr(ita_fisetin_data$PG.Genes, checkdbs)

pdf("C:/Users/piercetf/Documents/enrichr_fisetin_ita.pdf")

plotEnrich(fi[[1]], title="fisetin BP I")
plotEnrich(fi[[2]], title="fisetin CC I")
plotEnrich(fi[[3]], title="fisetin MF I")
plotEnrich(fi[[4]], title="fisetin R-ome I")
plotEnrich(fi[[5]], title="fisetin kegg I")
plotEnrich(fi[[6]], title="fisetin InterPro I")

dev.off()

nq <- enrichr(nparc_quercetin_data$PG.Genes, checkdbs)

pdf("C:/Users/piercetf/Documents/enrichr_quercetin_nparc.pdf")

plotEnrich(nq[[1]], title="quercetin BP N")
plotEnrich(nq[[2]], title="quercetin CC N")
plotEnrich(nq[[3]], title="quercetin MF N")
plotEnrich(nq[[4]], title="quercetin R-ome N")
plotEnrich(nq[[5]], title="quercetin kegg N")
plotEnrich(nq[[6]], title="quercetin InterPro N")

dev.off()

iq <- enrichr(ita_quercetin_data$PG.Genes, checkdbs)

pdf("C:/Users/piercetf/Documents/enrichr_quercetin_ita.pdf")

plotEnrich(iq[[1]], title="quercetin BP I")
plotEnrich(iq[[2]], title="quercetin CC I")
plotEnrich(iq[[3]], title="quercetin MF I")
plotEnrich(iq[[4]], title="quercetin R-ome I")
plotEnrich(iq[[5]], title="quercetin kegg I")
plotEnrich(iq[[6]], title="quercetin InterPro I")

dev.off()
