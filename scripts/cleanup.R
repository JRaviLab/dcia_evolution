library(tidyverse)
setwd("/data/research/jravilab/dcia_mtb/data/")
df <- read_tsv("DciA_combined_cln.tsv", col_names = T)
# get rid of headers from integration
df <- df[!(df$Query=="Query"),]
# Label DciA or DnaB
df <- df %>% mutate(AccNum = if_else(grepl("helicase", STitle), paste0(AccNum, "_DnaB"), paste0(AccNum, "_DciA")))
# seperate DNAB
df2 <- subset(df, grepl("helicase", STitle))
# seperate DciA/DUF721
df3 <- subset(df, ! grepl("helicase", STitle))
# Same for IPR
in_ipr <- read_tsv("DciA_combined_ipr.tsv", col_names = T)

in_ipr <- in_ipr[!(in_ipr$DB.ID=="DB.ID"),]
in_ipr <- in_ipr %>% mutate(AccNum = if_else(AccNum %in% df2$AccNum, paste0(AccNum, "_DnaB"), paste0(AccNum, "_DciA")))
ipr_2 <- subset(in_ipr, AccNum %in% df2$AccNum)
ipr_3 <- subset(in_ipr, AccNum %in% df3$AccNum)

write_tsv(df, "all_combined_cln.tsv")
write_tsv(df2, "DnaB_combined_cln.tsv")
write_tsv(df3, "DciA_combined_cln.tsv")

write_tsv(in_ipr, "all_combined_ipr.tsv")
write_tsv(ipr_2, "DnaB_combined_ipr.tsv")
write_tsv(ipr_3, "DciA_combined_ipr.tsv")