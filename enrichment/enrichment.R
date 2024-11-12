
library(org.Hs.eg.db)
library(clusterProfiler)

# get all candidate genes
#candidate_filename <- "T:/TGB/LSS/TGU/Users/Tomas/2024_CETSA_MS/Monocyte_CETSA_Statistical_Analysis/CETSA_ind_temp_analysis_starting_files/Candidates.tsv"
candidate_filename <- "C:/Users/piercetf/Downloads/CachedCETSAData/Candidates.tsv"
candidates <- read.csv(candidate_filename, sep="\t")
candidate_genes <- base::unique(candidates$Genes)
candidate_prots <- base::unique(candidates$UniProtIds)


# get nparc-detected significant genes for fisetin
fisetin_name <- "C:/Users/piercetf/Documents/nparc_unshared_fisetin_Oct2024.csv"
fisetin_table <- read.csv(fisetin_name)
fisetin_nparc_genes <- base::unique(fisetin_table$PG.Genes)
fisetin_nparc_prots <- base::unique(fisetin_table$PG.ProteinAccessions)

fisetin_nparc_overrep <- enrichGO(fisetin_nparc_prots, 
         "org.Hs.eg.db", 
         keyType="UNIPROT",
         ont="MF",
         universe=candidate_prots)

dotplot(fisetin_nparc_overrep)

fisetin_nparc_overrep <- enrichGO(fisetin_nparc_prots, 
                                  "org.Hs.eg.db", 
                                  keyType="UNIPROT",
                                  ont="CC",
                                  universe=candidate_prots)

dotplot(fisetin_nparc_overrep)

fisetin_nparc_overrep <- enrichGO(fisetin_nparc_prots, 
                                  "org.Hs.eg.db", 
                                  keyType="UNIPROT",
                                  ont="BP",
                                  universe=candidate_prots)

dotplot(fisetin_nparc_overrep)

# same idea but for quercetin
quercetin_name <- "C:/Users/piercetf/Documents/nparc_unshared_quercetin_Oct2024.csv"
quercetin_table <- read.csv(quercetin_name)
quercetin_nparc_prots <- base::unique(quercetin_table$PG.ProteinAccessions)

quercetin_nparc_overrep <- enrichGO(quercetin_nparc_prots, 
                                  "org.Hs.eg.db", 
                                  keyType="UNIPROT",
                                  ont="MF",
                                  universe=candidate_prots)

dotplot(quercetin_nparc_overrep)

quercetin_nparc_overrep <- enrichGO(quercetin_nparc_prots, 
                                    "org.Hs.eg.db", 
                                    keyType="UNIPROT",
                                    ont="CC",
                                    universe=candidate_prots)

dotplot(quercetin_nparc_overrep)

quercetin_nparc_overrep <- enrichGO(quercetin_nparc_prots, 
                                    "org.Hs.eg.db", 
                                    keyType="UNIPROT",
                                    ont="BP",
                                    universe=candidate_prots)

dotplot(quercetin_nparc_overrep)


# same approach but with ITA data
fisetin_name2 <- "C:/Users/piercetf/Documents/ITA_fisetin_Student_Oct2024.csv"

fisetin_table <- read.csv(fisetin_name2)
fisetin_nparc_genes <- base::unique(fisetin_table$PG.Genes)
fisetin_nparc_prots <- base::unique(fisetin_table$PG.ProteinAccessions)

fisetin_nparc_overrep <- enrichGO(fisetin_nparc_prots, 
                                  "org.Hs.eg.db", 
                                  keyType="UNIPROT",
                                  ont="MF",
                                  universe=candidate_prots)

dotplot(fisetin_nparc_overrep)

fisetin_nparc_overrep <- enrichGO(fisetin_nparc_prots, 
                                  "org.Hs.eg.db", 
                                  keyType="UNIPROT",
                                  ont="CC",
                                  universe=candidate_prots)

dotplot(fisetin_nparc_overrep)

fisetin_nparc_overrep <- enrichGO(fisetin_nparc_prots, 
                                  "org.Hs.eg.db", 
                                  keyType="UNIPROT",
                                  ont="BP",
                                  universe=candidate_prots)

dotplot(fisetin_nparc_overrep)

# and again for quercetin with individual temperatures
quercetin_name2 <- "C:/Users/piercetf/Documents/ITA_quercetin_Student_Oct2024.csv"
quercetin_table <- read.csv(quercetin_name2)
quercetin_nparc_prots <- base::unique(quercetin_table$PG.ProteinAccessions)

quercetin_nparc_overrep <- enrichGO(quercetin_nparc_prots, 
                                    "org.Hs.eg.db", 
                                    keyType="UNIPROT",
                                    ont="MF",
                                    universe=candidate_prots)

dotplot(quercetin_nparc_overrep)

quercetin_nparc_overrep <- enrichGO(quercetin_nparc_prots, 
                                    "org.Hs.eg.db", 
                                    keyType="UNIPROT",
                                    ont="CC",
                                    universe=candidate_prots)

dotplot(quercetin_nparc_overrep)

quercetin_nparc_overrep <- enrichGO(quercetin_nparc_prots, 
                                    "org.Hs.eg.db", 
                                    keyType="UNIPROT",
                                    ont="BP",
                                    universe=candidate_prots)

dotplot(quercetin_nparc_overrep)

