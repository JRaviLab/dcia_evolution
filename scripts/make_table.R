library(gt)
library(here)
library(tidyverse)
library(paletteer)
df <- read_tsv("/data/scratch/janani/molevolvr_out/RlSBhA_phylo/cln_combined.tsv")
df2 <- read_tsv("Helen_table.tsv")
df <- merge(df, df2, by="AccNum", all.x=TRUE)
df <- df %>% subset( select = c(Name, Species, TaxID, Lineage, DomArch.Pfam, DomArch.MobiDBLite,
Length,DUF721_range, Group, Gram_stain
))
df <- df %>% mutate(Species = if_else(Name == "BSpiroc_Linterrogans_WP_000650726.1", "Leptospira interrogans", Species))
df <- df %>% mutate(Species = if_else(Name == "BProteo_Babortus_WP_002963653.1", "Brucella abortus", Species))
df <- df %>% mutate(Species = if_else(Name == "BProteo_Paeruginosa_WP_003120896.1", "Pseudomonas aeruginosa", Species))
df <- df %>% mutate(Species = if_else(Name == "BActino_Scoelicolor_WP_003975057.1", "Streptomyces coelicolor", Species))
df <- df %>% mutate(Species = if_else(Name == "BProteo_Pmirabilis_WP_004244106.1", "Proteus mirabilis", Species))
df <- df %>% mutate(Species = if_else(Name == "BActino_Rjostii_WP_009476748.1", "Rhodococcus jostii", Species))
df <- df %>% mutate(Species = if_else(Name == "BActino_Savermiltilis_WP_010985745.1", "Streptomyces avermiltilis", Species))
df <- df %>% mutate(Species = if_else(Name == "BProteo_Maustralicum_WP_015318768.1", "Mesorhizobium australicum", Species))
df <- df[!duplicated(df$Name), ]
table <- df %>%
  gt() %>%
  cols_label(Gram_stain = "Gram stain", DUF721_range = "DUF721 range (aa)", Group = "DciA group") %>%
  tab_header(
    title=md("Table 1: DciA query proteins used to identify homologs across the bacterial kingdoms")) %>%
  tab_options(
      # Headings; Titles
      heading.background.color="black",
      heading.border.bottom.color="#989898",
      heading.title.font.size="12px",
      heading.subtitle.font.size="11px",
      # Column labels
      column_labels.background.color="grey50", #B09C85FF
      column_labels.font.size="12px",
      # Stubs
      stub.background.color="#4DBBD5", #B09C85FF
      stub.border.style="dashed",
      stub.border.color="#989898",
      stub.border.width="1px",
      # Row groups
      row_group.background.color="#3C5488", #FFEFDB80
      row_group.border.top.color="#989898",
      row_group.border.bottom.style="none",
      row_group.font.size="12px",
      # Summary rows
      summary_row.border.color="#989898",
      # summary_row.background.color="#FFEBEE",
      # grand_summary_row.background.color="#FFFFFF",
      # Table
      table.font.color="#323232",
      table_body.hlines.color="#989898",
      table_body.border.top.color="#989898",
      table.font.size="10px",
      table.width="80%"
    )

gtsave(table, "table.pdf")
