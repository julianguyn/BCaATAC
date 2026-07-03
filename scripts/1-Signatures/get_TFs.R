# helper script to grab TFs for chromvar

suppressPackageStartupMessages({
    library(tidyverse)
    library(readr)
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
})

source("utils/palettes.R")

set.seed(101)

###########################################################
# Load in HOMER results
###########################################################

homer <- data.frame(matrix(nrow=0, ncol=13))

for (arche in paste0("ARCHE", 1:6)) {
  filename <- paste0("data/results/data/1-Signatures/findMotifsGenome/", arche, "/knownResults.txt")
  df <- read_tsv(filename, show_col_types = FALSE)
  colnames(df) <- c(
    "motif_name", "consensus", "pvalue", "log_pvalue", "qvalue",
    "num_target_seqs", "pct_target_seqs",
    "num_bg_seqs", "pct_bg_seqs"
  )

  df <- df %>%
    mutate(
      pct_target_seqs = as.numeric(gsub("%", "", pct_target_seqs)),
      pct_bg_seqs     = as.numeric(gsub("%", "", pct_bg_seqs)),
      pct_diff        = pct_target_seqs - pct_bg_seqs,
      fold_enrichment = pct_target_seqs / pct_bg_seqs
    )

  df <- df %>%
    mutate(tf_name = sub("\\(.*", "", motif_name))

  df$ARCHE <- arche
  #df <- df[df$qvalue < 0.05,]
  homer <- rbind(homer, df) |> as.data.frame()
}

# filter TFs
filtered <- homer %>%
    filter(
      qvalue < 0.01,
      pct_target_seqs > 10,      # >10% presence
      fold_enrichment > 1.5     # >1.25x background rate
    ) %>%
    arrange(desc(pct_diff))
homer <- homer[homer$motif_name %in% filtered$motif_name,]

# label the FOXA1 dups
homer$tf_name[grepl("LNCAP-FOXA1-ChIP-Seq", homer$motif_name, fixed = TRUE)] <- "FOXA1 (LNCAP)"
homer$tf_name[grepl("MCF7-FOXA1-ChIP-Seq", homer$motif_name, fixed = TRUE)]  <- "FOXA1 (MCF7)"

###########################################################
# Plot dataframe
###########################################################

toPlot <- homer %>% 
  select(ARCHE, tf_name, fold_enrichment) %>%
  pivot_wider(
    names_from = ARCHE,
    values_from = fold_enrichment
  ) %>% 
  column_to_rownames("tf_name") %>% as.data.frame()

# create annotation table
tf_annotation <- data.frame(
  tf_name = rownames(toPlot),
  family = c(
    "Zf (CTCF)", "Sp/KLF", "ETS", "ETS", "ETS",
    "ETS", "Zf (CTCF)", "bZIP (AP-1)", "bZIP (AP-1)", "bZIP (AP-1)",
    "bZIP (AP-1)", "bZIP (AP-1)", "bZIP (AP-1)", "bZIP (AP-1)", "CBF/NF-Y",
    "bZIP (AP-1)", "bZIP (AP-1)", "ETS", "Grainyhead", "Forkhead",
    "ETS", "IRF", "ETS", "Forkhead", "IRF"
  ),
  immune_lineage = c(
    FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, TRUE,  FALSE, FALSE,
    FALSE, TRUE,  TRUE,  FALSE, TRUE
  ),
  stringsAsFactors = FALSE
)

###########################################################
# Set colours
###########################################################

family_levels <- unique(tf_annotation$family)
family_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(family_levels)),
  family_levels
)

immune_colors <- c("TRUE" = "#B56576", "FALSE" = "grey85")

cols <- colorRampPalette(c("#466D9F", "#6878C9", "#7986C7", "#8DC2A4", "#EBE4C5", "#EAEAEA"))(9)
col_fun <- colorRamp2(
    seq(3,0,
        length.out = 9),
    cols
)

###########################################################
# Create annotations
###########################################################

col_ha <- HeatmapAnnotation(
    'ARCHE' = paste0("ARCHE", 1:6),
    col = list('ARCHE' = ARCHE_pal),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

row_ha <- rowAnnotation(
  Family = tf_annotation$family,
  `Immune Lineage` = tf_annotation$immune_lineage,
  col = list(
    Family = family_colors,
    `Immune Lineage` = immune_colors
  ),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    Family = list(title = "TF Family", title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    `Immune Lineage` = list(title = "Immune Lineage TF", title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8))
  )
)

###########################################################
# Heatmap
###########################################################

ht <- Heatmap(
  toPlot,
  row_split = 6,
  cluster_columns = FALSE,
  name = "Fold\nEnrichment",
  col = col_fun,
  row_names_gp = gpar(fontsize = 8),
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 0,
  column_names_centered = TRUE,
  row_title = "Transcription Factor",
  row_title_side = "left",
  row_title_rot = 90,
  row_title_gp = gpar(fontsize = 9),
  bottom_annotation = col_ha,
  right_annotation = row_ha
)
filename <- "data/results/figures/1-Signatures/TF_heatmap.png"
png(filename, width = 6.5, height = 5, res = 600, units = "in")
ht
dev.off()




###########################################################
# OLD CODE
###########################################################


library(readxl)

# make dataframe to store results
tf <- data.frame(matrix(nrow=0, ncol=7))

for (arche in paste0("ARCHE", 1:6)) {
    label <- paste(gsub("RCHE", "", arche), "K", sep = " - ")

    file <- "data/results/data/1-Signatures/findMotifsGenome/findMotifsGenome.xlsx"
    motifs <- read_excel(file, sheet = label) |> as.data.frame()

    # get variables and format
    motifs$Name <- sub("\\/.*", "", motifs$Name)
    colnames(motifs) <- gsub(" Sequences with Motif", "", colnames(motifs))
    motifs <- motifs[,c("Name", "Rank", "P-value", "% of Targets", "% of Background")]
    motifs$Rank <- factor(motifs$Rank, levels = motifs$Rank)

    # compute fold change
    motifs$FoldChange <- motifs$'% of Targets' / motifs$'% of Background'

    # split name
    motifs$family <- sub("\\)", "", sub(".*\\(", "", motifs$Name))
    motifs$Name <- sub("\\(.*", "", motifs$Name)

    motifs$ARCHE <- arche

    tf <- rbind(tf, motifs)
}

# get unique genes
genes <- unique(tf$Name)

motif_to_gene <- list(
  "Sp1"           = "SP1",
  "KLF17"         = "KLF17",
  "ETS"           = NA,             # family motif, not a single gene
  "Elk4"          = "ELK4",
  "KLF3"          = "KLF3",
  "Klf15"         = "KLF15",
  "NFY"           = NA,             # complex of NFYA/NFYB/NFYC, not a single gene
  "Elk1"          = "ELK1",
  "Sp5"           = "SP5",
  "ELF1"          = "ELF1",
  "KLF1"          = "KLF1",
  "Fli1"          = "FLI1",
  "GFY-Staf"      = "ZNF143",       # composite motif; Staf component = ZNF143, GFY unresolved
  "GABPA"         = "GABPA",
  "GRHL2"         = "GRHL2",
  "NF1-halfsite"  = "NFIX",         # half-site variant of NF1 family motif; ambiguous among NFIA/NFIB/NFIC/NFIX
  "Sox3"          = "SOX3",
  "Sox9"          = "SOX9",
  "Sox2"          = "SOX2",
  "Sox10"         = "SOX10",
  "Sox17"         = "SOX17",
  "Sox15"         = "SOX15",
  "NF1"           = NA,             # family motif; ambiguous among NFIA/NFIB/NFIC/NFIX
  "Sox4"          = "SOX4",
  "Sox21"         = "SOX21",
  "Sox6"          = "SOX6",
  "SOX1"          = "SOX1",
  "Sox7"          = "SOX7",
  "Tlx?"          = NA,             # HOMER flags as low-confidence; could be TLX1/TLX2/NR2E1
  "FOXM1"         = "FOXM1",
  "FOXA1"         = "FOXA1",
  "Foxa2"         = "FOXA2",
  "Fox:Ebox"      = NA,             # composite motif (FOX site + E-box), not a single gene
  "Foxa3"         = "FOXA3",
  "AP-2alpha"     = "TFAP2A",
  "AP-2gamma"     = "TFAP2C",
  "FoxL2"         = "FOXL2",
  "PHA-4"         = "FOXA2",        # C. elegans gene; closest mammalian ortholog (FOXA family)
  "FoxD3"         = "FOXD3",
  "Foxo3"         = "FOXO3",
  "FOXK1"         = "FOXK1",
  "Foxf1"         = "FOXF1",
  "FOXP1"         = "FOXP1",
  "FOXK2"         = "FOXK2",
  "CTCF"          = "CTCF",
  "PU.1"          = "SPI1",
  "ETS1"          = "ETS1",
  "SpiB"          = "SPIB",
  "Elf4"          = "ELF4",
  "Etv2"          = "ETV2",
  "ETV1"          = "ETV1",
  "ERG"           = "ERG",
  "Ets1-distal"   = "ETS1",         # distal genomic site variant of the ETS1 motif
  "ELF5"          = "ELF5",
  "IRF8"          = "IRF8",
  "ETV4"          = "ETV4",
  "ELF3"          = "ELF3",
  "BORIS"         = "CTCFL"
)

unique_genes <- unique(na.omit(unlist(motif_to_gene)))
writeLines(unique_genes, "data/procdata/TCGA/homer_TF_genes.txt")

motifs_df <- data.frame(
  motif = names(motif_to_gene),
  gene = unlist(motif_to_gene)
)
tf$gene <- motifs_df$gene[match(tf$Name, motifs_df$motif)]

write.table(tf, "data/procdata/TCGA/homer_compiled_TF_genes.txt")