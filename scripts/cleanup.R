library(tidyverse)
setwd("/data/research/jravilab/dcia_mtb/data/")
df <- read_tsv("DciA_combined_cln.tsv", col_names = T)
# get rid of headers from integration
df <- df[!(df$Query=="Query"),]
# Label DciA or DnaB
df <- df %>% mutate(QueryName = if_else(grepl("helicase", STitle), paste0(QueryName, "_DnaB"), paste0(QueryName, "_DciA")))
# seperate DNAB
df2 <- subset(df, grepl("helicase", STitle))
# seperate DciA/DUF721
df3 <- subset(df, ! grepl("helicase", STitle))
# Same for IPR
in_ipr <- read_tsv("DciA_combined_ipr.tsv", col_names = T)

in_ipr <- in_ipr[!(in_ipr$DB.ID=="DB.ID"),]
#in_ipr <- in_ipr %>% mutate(QueryName = if_else(AccNum %in% df2$AccNum, paste0(QueryName, "_DnaB"), paste0(QueryName, "_DciA")))
ipr_2 <- subset(in_ipr, AccNum %in% df2$AccNum)
ipr_3 <- subset(in_ipr, AccNum %in% df3$AccNum)

write_tsv(df, "/data/scratch/janani/molevolvr_out/Dciacb_full/cln_combined.tsv")
write_tsv(df2, "/data/scratch/janani/molevolvr_out/DnaBcm_full/cln_combined.tsv")
write_tsv(df3, "/data/scratch/janani/molevolvr_out/DciAcm_full/cln_combined.tsv")

write_tsv(in_ipr, "/data/scratch/janani/molevolvr_out/Dciacb_full/ipr_combined.tsv")
write_tsv(ipr_2, "/data/scratch/janani/molevolvr_out/DnaBcm_full/ipr_combined.tsv")
write_tsv(ipr_3, "/data/scratch/janani/molevolvr_out/DciAcm_full/ipr_combined.tsv")