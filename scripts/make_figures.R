## Installing dependencies
library(tidyverse)
library(aplot)
library(ape)
#library(ggtree)
library(tidytree)
library(seqinr)
library(msa)
library(htmltools)
library(data.table)


## Sourcing relevant scripts from MolEvolvR back-end
## Burke*, Chen*, Sosinski*, et al., bioRxiv 2022 | http://jravilab.org/molevolvr
source("/data/research/jravilab/molevol_scripts/R/cleanup.R")
source("/data/research/jravilab/molevol_scripts/R/summarize.R")
source("/data/research/jravilab/molevol_scripts/R/plotting.R")
source("/data/research/jravilab/molevol_scripts/R/networks_domarch.R")
source("/data/research/jravilab/molevol_scripts/R/networks_gencontext.R")
source("/data/research/jravilab/molevol_scripts/R/pre-msa-tree.R")
source("/data/research/jravilab/molevol_scripts/R/ipr2viz.R")
source("/data/research/jravilab/molevol_scripts/R/lineage.R")
source("/data/research/jravilab/molevol_scripts/R/msa.R")
#source("/data/research/jravilab/molevol_scripts/scripts/tree.R")
source("/data/research/jravilab/molevol_scripts/R/colnames_molevol.R")
source("/data/research/jravilab/molevolvr_app/pins/PinGeneration.R")
source("/data/research/jravilab/molevolvr_app/scripts/ui/components.R")
source("/data/research/jravilab/molevolvr_app/scripts/ui/UIOutputComponents.R")
source("/data/research/jravilab/molevolvr_app/scripts/ui/splashPageComponent.R")
source("/data/research/jravilab/molevolvr_app/scripts/MolEvolData_class.R")
source("/data/research/jravilab/molevolvr_app/scripts/ui/tabText.R")
source("/data/research/jravilab/molevolvr_app/scripts/utility.R")
source("/data/research/jravilab/molevol_data/scripts/dcia_figures.R")
conflicted::conflict_prefer("count", "dplyr")

## Importing data
df <- read_tsv("/data/scratch/jburke/dcia/all_seqs_full_with_mobi_ipr_distinct_2.tsv")
df$Name <- str_replace_all(df$Name, "\\[", "")
df$QueryName <- df$Name
in_ipr <- read_tsv("/data/scratch/jburke/dcia/all_seqs.iprscan_cln.tsv")

## FIGURE 1
#df for Fig 1 filter out those with no lineage
df <- df %>% filter(Lineage != "") %>% drop_na(Lineage)
df <- df %>% mutate(DomArch.Pfam = gsub("_[0-9]{1,2}", "", DomArch.Pfam)) %>% mutate(DomArch.Pfam = str_replace(DomArch.Pfam, "TPR\\+.*", paste0("TPR(", str_count(DomArch.Pfam, "TPR"), ")"))) %>%
  mutate(DomArch.Pfam = str_replace(DomArch.Pfam, "DciA\\+DciA\\+DciA\\+DciA", "DciA(4)"))

# Fig 1a -> sunburst of everything without metagenomes -> 1 bact 1 non bact
df_all <- df %>% filter(!grepl("meta", Organism))
df_bact <- df %>% filter(grepl("Bacteria>", Lineage))
df_other <- df %>% filter(!grepl("Bacteria>", Lineage))

# Fig 1b -> Heatmap of everything without metagenomes 0> maybe 1 bact/non bact -> remove dcia alone
# sunburst code
Lineage_domarch <- lineage.DA.plot2(df_all, colname = "DomArch.Pfam", cutoff = 100, RowsCutoff = F, text_size = 12)
ggsave("dcia_figs/fig1/1B_all.png", Lineage_domarch, dpi = 400, device = "png", height = 10, width = 25)
df_dcia <- df_all %>%
  filter(DomArch.Pfam == "DciA" & grepl(">", Lineage))
df_dcia$Lineage <- gsub("Bacteria>Candid.*", "Bacteria>Candidatus", df_dcia$Lineage)
df_dcia$Lineage <-  gsub("Bacteria>candid.*", "Bacteria>Candidatus", df_dcia$Lineage)
df_dcia$Lineage <-  gsub("Archaea>Candid.*", "Archaea>Candidatus", df_dcia$Lineage)

Lineage_domarch <- lineage.DA.plot2(df_dcia, colname = "DomArch.Pfam", cutoff = 100, RowsCutoff = F, text_size = 14)
ggsave("dcia_figs/fig1/1B_dcia.png", Lineage_domarch, dpi = 400, device = "png", height = 4, width = 12)

df_no_dcia <- df_all %>%
  filter(DomArch.Pfam != "DciA") %>%
  filter(grepl("Bacteria>", Lineage))
Lineage_domarch <- lineage.DA.plot2(df_no_dcia, colname = "DomArch.Pfam",
                                    cutoff = 100, RowsCutoff = F, text_size = 14)
ggsave("dcia_figs/fig1/1B_no_dcia.png", Lineage_domarch,
       dpi = 400, device = "png", height = 8, width = 8)

# Fig 1C -> unique domarch of only bacteria + tree + MSA
df_bact <- df_bact %>%
  filter(grepl("Bacteria", Lineage)) %>%
  drop_na(AccNum) %>%
  filter(!grepl("candidate", Lineage)) %>%
  filter(!grepl("Candidatus", Lineage)) %>%
  filter(!grepl("candidate", Species)) %>%
  filter(!grepl("Candidat", Species)) %>%
  filter(!grepl("uncultur", Species))

df_bact$DomArch.Pfam <- str_replace_all(df_bact$DomArch.Pfam , "Hydrolase\\+DciA", "HAD_2+DciA")
df_unique <- df_bact %>% distinct(DomArch.Pfam, .keep_all = TRUE)
df_unique <- df_unique %>% filter(Accession != "A0A2G4H134")
in_ipr_bact <- in_ipr
in_ipr_bact$SignDesc <- str_replace_all(in_ipr_bact$SignDesc,
                                        "consensus disorder prediction", "Intrinsically disordered region (MobiDBLite)")
in_ipr_bact$SignDesc <- str_replace_all(in_ipr_bact$SignDesc,
                                        "haloacid dehalogenase-like hydrolase", "Haloacid dehalogenase-like hydrolase")
in_ipr_bact$SignDesc <- str_replace_all(in_ipr_bact$SignDesc,
                                        "Tetratricopeptide repeat", "TPR repeat")
in_ipr_bact$ShortName <- str_replace_all(in_ipr_bact$ShortName,
                                         "Hydrolase", "HAD_2")

domarch <- ipr2viz_web2(infile_ipr = in_ipr_bact, accessions = df_unique$Name,
                        analysis = c("Domarch"), group_by = "Analysis", text_size = 16)
ggsave("dcia_figs/fig1/1C_slim.png",
       domarch, dpi = 400, device = "png", height = 12, width = 10)
ggsave("dcia_figs/fig1/1C_full.png",
       domarch, dpi = 400, device = "png", height = 12, width = 20)

## FIGURE S1
# Fig S1a Domarch of 2 euk reps
df_euks <- df %>% filter(AccNum == "KAF6150485.1" | AccNum == "PMD57303.1")
write_tsv(df_euks, "dcia_figs/figS1/table.tsv")
domarch <- ipr2viz_web2(infile_ipr = in_ipr, accessions = df_euks$Name,
                        analysis = c("Pfam"), group_by = "Analysis", text_size = 14, cols = 2, rows = 1)
ggsave("dcia_figs/figS1/S1A.png", domarch, dpi = 400, device = "png", height = 4, width = 8)

## Table S2
# Table S1 | archaea and eukaryotes
df_others <- df %>% filter(!grepl("Bacteria", Lineage))
write_tsv(df_others, "dcia_figs/T1/table.tsv")

## Supplementary FIGURE S2
# Fig S2a -> sunburst of starting points
df <- read_tsv("/data/scratch/janani/molevolvr_out/v134LF_full/cln_combined_original.tsv")
ipr <- read_tsv("/data/scratch/janani/molevolvr_out/v134LF_full/ipr_combined.tsv")
df <- df %>% filter(grepl("DciA", DomArch.Pfam)) %>%
  filter(grepl("Bacteria>", Lineage)) %>%
  filter(!grepl("Candidatus", Lineage)) %>%
  arrange(desc(PcPositive)) %>%
  distinct(AccNum, .keep_all = TRUE)
ipr <- ipr %>% filter(ipr$AccNum %in% df$AccNum)
fwrite(df, "/data/scratch/janani/molevolvr_out/v134LF_full/cln_combined.tsv",
       sep = "\t", quote = FALSE)
fwrite(ipr,"/data/scratch/janani/molevolvr_out/v134LF_full/ipr_combined.tsv",
       sep = "\t", quote = FALSE)

# [REMOVED] Fig S2x Domarch heatmap of molevolvr homologs
df <- df %>%
  mutate(DomArch.Pfam = gsub("_[0-9]{1,2}", "", DomArch.Pfam)) %>%
  mutate(DomArch.Pfam = str_replace(DomArch.Pfam, "TPR\\+.*",
                                    paste0("TPR(", str_count(DomArch.Pfam, "TPR"), ")"))) %>%
  mutate(DomArch.Pfam = str_replace(DomArch.Pfam, "DciA\\+DciA", "DciA(2)"))%>%
  mutate(DomArch.Pfam = str_replace(DomArch.Pfam, "DUF2726\\+DciA", "TPR(1)+DUF2726+DciA"))

Lineage_domarch <- lineage.DA.plot2(df, colname = "DomArch.Pfam", cutoff = 100, RowsCutoff = F, text_size = 16)
ggsave("dcia_figs/figS2/S2B.png", Lineage_domarch,
       dpi = 400, device = "png", height = 7, width = 13)

# Fig S2b Domarch fig of unique domarchs from molevolvr
df <- df %>% mutate(DomArch.Pfam = str_replace_all(DomArch.Pfam, "Hydrolase", "HAD_2"))
df_unique <- df %>% distinct(DomArch.Pfam, .keep_all = TRUE)
ipr$SignDesc <- str_replace_all(ipr$SignDesc, "haloacid dehalogenase-like hydrolase", "Haloacid dehalogenase-like hydrolase")
ipr$SignDesc <- str_replace_all(ipr$SignDesc, "Tetratricopeptide repeat", "TPR repeat")
ipr$SignDesc <- str_replace_all(ipr$SignDesc, "consensus disorder prediction", "Intrinsic disorder region (MobiDBLite)")
ipr$ShortName <- str_replace_all(ipr$ShortName, "consensus disorder prediction", "Intrinsic disorder region (MobiDBLite)")

ipr_plot <- ipr2viz_web2(infile_ipr = ipr, accessions = df_unique$Name,
                         analysis = c("Domarch"), group_by = "Analysis",
                         rows = 8, cols = 2, text_size = 16)
ggsave("dcia_figs/figS2/S2C_full.png", ipr_plot,
       dpi = 400, device = "png", height = 6, width = 14)
ipr_plot <- ipr2viz_web2(infile_ipr = ipr, accessions = df_unique$Name,
                         analysis = c("Domarch"), group_by = "Analysis",
                         rows = 8, cols = 2, text_size = 16)
ggsave("dcia_figs/figS2/S2C_slim.png", ipr_plot,
       dpi = 400, device = "png", height = 6, width = 10)

## Figure 2
# Fig 2a length box plot of DciA only bac
df <- read_tsv("/data/scratch/jburke/dcia/all_seqs_full_with_mobi_ipr_distinct_2.tsv")
df$Name <- str_replace_all(df$Name, "\\[", "")
df$QueryName <- df$Name
in_ipr <- read_tsv("/data/scratch/jburke/dcia/all_seqs.iprscan_cln.tsv")
# df for fig 1 filter out those with no lineage
df <- df %>% filter(Lineage != "") %>% drop_na(Lineage)
df <- df %>% mutate(DomArch.Pfam = gsub("_[0-9]{1,2}", "", DomArch.Pfam)) %>%
  mutate(DomArch.Pfam = str_replace(DomArch.Pfam, "TPR\\+.*",
                                    paste0("TPR(", str_count(DomArch.Pfam, "TPR"), ")"))) %>%
  mutate(DomArch.Pfam = str_replace(DomArch.Pfam, "DciA\\+DciA\\+DciA\\+DciA", "DciA(4)"))
df_fig2 <- df %>%
  filter(DomArch.Pfam == "DciA") %>%
  filter(grepl("Bacteria>", Lineage)) %>%
  filter(!grepl("Candidatus", Lineage)) %>%
  filter(!grepl("candidate", Lineage)) %>%
  distinct(AccNum, .keep_all = TRUE) %>%
  drop_na(pfam_start)

df_fig2$domainLength <- df_fig2$pfam_stop - df_fig2$pfam_start
df_fig2$fromStart <- df_fig2$pfam_start
df_fig2$fromEnd <- df_fig2$seq_length - df_fig2$pfam_stop
df_fig2$seqLength <- df_fig2$seq_length
df_fig2 <- df_fig2 %>%
  pivot_longer(c(seqLength, fromStart, fromEnd), names_to = "Type", values_to = "length")
df_fig2$Type <- factor(df_fig2$Type, levels = c("seqLength", "fromStart", "fromEnd"))
df_fig2_stats <- df_fig2 %>%
  group_by(Type) %>%
  mutate(med = median(length), twent = quantile(length, c(0.25)), sev = quantile(length, c(0.75)))
fig_2a <- ggplot(df_fig2, aes(x=Lineage, length)) +
  geom_boxplot(outlier.alpha = 0.5, outlier.color = "black") +
  labs(y = "Length", x = "Lineage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.background = element_rect("white", size = 0),
        text = element_text(size = 18), panel.border = element_blank()) +
  facet_grid(Type ~ . , scales = "free") +
  geom_hline(data = df_fig2_stats, aes(yintercept = med), color = "blue") +
  geom_hline(data = df_fig2_stats, aes(yintercept = twent), color = "blue", alpha = 0.5) +
  geom_hline(data = df_fig2_stats, aes(yintercept = df_fig2_stats$sev), color = "blue", alpha = 0.5)
ggsave("dcia_figs/fig2/2a.png", fig_2a,
       dpi = 400, device = "png", height = 14, width = 16)

## FIGURE 3
# Fig 3a Mobidb heatmap on only bacteria
Lineage_domarch <- lineage.DA.plot2(df_bact, colname = "DomArch.MobiDBLite",
                                    cutoff = 100, RowsCutoff = F, text_size = 12)
ggsave("dcia_figs/fig3/3A.png", Lineage_domarch,
       dpi = 400, device = "png", height = 4, width = 8)

df_fig3 <- df %>% filter(grepl("DciA", DomArch.Pfam)) %>%
  filter(grepl("consensus disorder prediction", DomArch.MobiDBLite)) %>%
  filter(grepl("Bacteria>", Lineage)) %>%
  drop_na(AccNum) %>%
  filter(!grepl("candidate", Lineage)) %>%
  filter(!grepl("Candidatus", Lineage)) %>%
  filter(!grepl("candidate", Species)) %>%
  filter(!grepl("Candidat", Species)) %>%
  filter(!grepl("uncultur", Species)) %>%
  drop_na(mobi_stop)

df_fig3$domainLength <- df_fig3$mobi_stop - df_fig3$mobi_start
df_fig3$fromStart <- df_fig3$mobi_start
df_fig3$fromEnd <- df_fig3$seq_length - df_fig3$mobi_stop
df_fig3 <- df_fig3 %>%
  pivot_longer(c(domainLength, fromStart, fromEnd), names_to = "Type", values_to = "length")
df_fig3$Type <- factor(df_fig3$Type, levels = c("domainLength", "fromStart", "fromEnd"))
df_fig3_stats <- df_fig3 %>% group_by(Type) %>%
  mutate(med = median(length), twent = quantile(length, c(0.25)), sev = quantile(length, c(0.75)))

b <- ggplot(df_fig3, aes(x=Lineage, length)) +
  geom_boxplot(outlier.alpha = 0.5, outlier.color = "black") +
  labs(y = "Length", x = "Lineage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.background = element_rect("white", size = 0) ,text = element_text(size = 18),
        panel.border = element_blank()) +
  facet_grid(Type ~ . , scales = "free") +
  geom_hline(data = df_fig3_stats, aes(yintercept = med), color = "blue") +
  geom_hline(data = df_fig3_stats, aes(yintercept = twent), color = "blue", alpha = 0.5) +
  geom_hline(data = df_fig3_stats, aes(yintercept = sev), color = "blue", alpha = 0.5)
ggsave("dcia_figs/fig3/3_facet.png", b,
       dpi = 400, device = "png", height = 14, width = 16)

# Fig 3e group by mobidb, lineage, select representative, run MSA
df_fig3 <- df_fig3 %>%
  group_by(DomArch.MobiDBLite, Lineage) %>%
  distinct(Accession, .keep_all = TRUE)
write_tsv(df_fig3, "dcia_figs/fig3/reps.tsv")

## FIGURE 4
# Fig 4a stacked bar plot of lineage
df_fig4 <- df_bact %>% filter(DomArch.Pfam == "DciA")%>% drop_na(seq_length)
df_fig4$fromEnd <- df_fig4$seq_length - df_fig4$pfam_stop
df_fig4$fromStart <- df_fig4$pfam_start
df_fig4$group <- ""
df_fig4$group[which(df_fig4$fromStart >= 14 & df_fig4$fromEnd < 14)] <- 1
df_fig4$group[which(df_fig4$fromEnd >= 14 &  df_fig4$fromStart < 14)] <- 2
df_fig4$group[which(df_fig4$fromEnd >= 14 &  df_fig4$fromStart >= 14)] <- 3
df_fig4$group[which(df_fig4$group == "")] <- 4
write_tsv(df_fig4, "dcia_figs/fig4/cln_with_groups.tsv")
stacked <- stacked_lin_plot(df_fig4, column = "group", cutoff = 100,
                            xlabel = "Group", legend.position = c(0.7,0.3),
                            legend.cols = 2)
ggsave("dcia_figs/fig4/4a.png", stacked,
       dpi = 400, device = "png", height = 15, width = 15)
stacked <- stacked_lin_plot(df_fig4, column = "group", cutoff = 100,
                            xlabel = "Group", legend.position = c(0.7,0.3),
                            legend.cols = 2, legend = FALSE)
ggsave("dcia_figs/fig4/4_no_legend.png", stacked,
       dpi = 400, device = "png", height = 15, width = 15)

# Fig 4b
df_fig4_b <- df_fig4 %>% filter(Lineage == "Bacteria>Proteobacteria" | Lineage == "Bacteria>Actinobacteria" | Lineage == "Bacteria>Bacteroidetes")
stacked <- stacked_lin_plot(df_fig4_b, column = "group",
                            Lineage_col = "Lineage_short", cutoff = 100,
                            xlabel = "Group", legend.position = c(0.7,0.4),
                            legend.cols = 2, reduce_lineage = FALSE)
ggsave("dcia_figs/fig4/4b.png", stacked, dpi = 400, device = "png", height = 15, width = 15)
stacked <- stacked_lin_plot(df_fig4_b, column = "group",
                            Lineage_col = "Lineage_short", cutoff = 100,
                            xlabel = "Group", legend.position = c(0.7,0.4),
                            legend.cols = 2, reduce_lineage = FALSE, legend = FALSE)
ggsave("dcia_figs/fig4/4b_no_legend.png", stacked,
       dpi = 400, device = "png", height = 15, width = 15)

# Fig 4c
df_fig4_c <- df_fig4 %>%
  filter(Lineage != "Bacteria>Proteobacteria" & Lineage != "Bacteria>Actinobacteria" & Lineage != "Bacteria>Bacteroidetes")
stacked <- stacked_lin_plot(df_fig4_c, column = "group", cutoff = 100,
                            xlabel = "Group", legend.position = c(0.7,0.15),
                            legend.cols = 4)
ggsave("dcia_figs/fig4/4c.png", stacked,
       dpi = 400, device = "png", height = 15, width = 15)
stacked <- stacked_lin_plot(df_fig4_b, column = "group", cutoff = 100,
                            xlabel = "Group", legend.position = c(0.7,0.15),
                            legend.cols = 4, legend = FALSE)
ggsave("dcia_figs/fig4/4c_no_legend.png", stacked,
       dpi = 400, device = "png", height = 15, width = 15)
