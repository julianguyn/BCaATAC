# score REFLECT 6b samples for TF binding sites of interest

# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(readxl)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
  library(tidyverse)
  library(umap)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
})

source("utils/score_arche_cfDNA.R")
source("utils/palettes.R")

set.seed(101)

# -----------------------------------------------------------
# Parse args
# -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
analysis <- args[1]

valid <- c("NFE2LC", "MYC", "ATF2", "UGT1A6")
if (is.na(analysis) || !analysis %in% valid) {
  stop(
    sprintf("Invalid analysis argument '%s'. Must be one of: %s",
            analysis, paste(valid, collapse = ", ")),
    call. = FALSE
  )
}
dir <- paste0("TF_Scores_Mitchell/", analysis)

###########################################################
# Load in data
###########################################################

# load in metadata
meta <- read_excel("data/rawdata/cfDNA/REFLECT-unfiltered/REFLECT-6B_metadata.xlsx", sheet = 1) |> as.data.frame()

###########################################################
# Format metadata
###########################################################

meta$Subtype_final[meta$Subtype_final == "ER+/HER2-"] <- "ER+"
meta$Subtype_final[meta$Subtype_final == "TBNC"] <- "TNBC"
meta$id_6b <- gsub("_", "-", gsub("_LB.*", "", meta$library_id))
meta$Subtype_final[meta$Subtype_final == "HER2+"] <- "HER2"
meta$REF052 <- ifelse(meta$sample_id == "REF052", "REF052", "Other")

###########################################################
# Get TF scores
###########################################################

scores <- score_arche_cfDNA(dir, noARCHE = TRUE)

###########################################################
# Plot heatmap
###########################################################

# helper function
plot_heatmap <- function(scores, analysis) {

  df <- scores %>%
    pivot_wider(
      names_from = Label,
      values_from = Score
    ) %>%
    column_to_rownames("Sample") %>%
    t() %>% as.data.frame()

  # set colours for plotting
  lim <- max(c(abs(min(df)), max(df)))
  score_pal = colorRamp2(seq(-lim, lim, length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))
  
  ha <- HeatmapAnnotation(
      Subtype = meta$Subtype_final[match(colnames(df), meta$id_6b)],
      REF052 = meta$REF052[match(colnames(df), meta$id_6b)],
      col = list(Subtype = cfDNA_subtype_pal, REF052 = REF052_pal))

  filename <- paste0("data/results/figures/Misc/TF_Scoring/", analysis, "_heatmap.png")
  png(filename, width = 9, height = 3, res = 600, units = "in")
  print(
      Heatmap(df, cluster_rows = FALSE, name = "TF Binding\nScore", col = score_pal,
          column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
          row_names_gp = gpar(fontsize = 10), top_annotation = ha)
  )
  dev.off()
}

plot_heatmap(scores, analysis)

###########################################################
# Plot barplots
###########################################################

scores$REF052 <- ifelse(scores$Sample == "REFLECT-0052-03", "REF052", "Other")
scores <- scores[order(scores$Score, decreasing = TRUE),]
scores$Sample <- factor(scores$Sample, levels = unique(scores$Sample))

p <- ggplot(scores, aes(x = Sample, y = Score, fill = REF052)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = REF052_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm')
    )

if (analysis == "UGT1A6") {
  p <- ggplot(scores, aes(x = Sample, y = Score, fill = Label)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge()) +
    scale_fill_manual(values = binary_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm')
    ) + geom_text(
      data = scores %>% filter(Sample == "REFLECT-0052-03"),
      aes(label = "*"),
      position = position_dodge(width = 0.9),
      vjust = -0.7,
      size = 6
  )
}

filename <- paste0("data/results/figures/Misc/TF_Scoring/", analysis, "_scores.png")
png(filename, width=8, height=5, units='in', res = 600, pointsize=80)
p
dev.off()

###########################################################
# Plot coverage plots
# modify from utils/plots/cfdna.R
###########################################################

plot_coverage <- function(cov, meta, analysis) {
    toPlot <- reshape2::melt(cov)
    toPlot$subtype <- meta$Subtype_final[match(toPlot$sample, meta$id_6b)]
    toPlot$REF052 <- meta$REF052[match(toPlot$sample, meta$id_6b)]
    toPlot$variable <- as.numeric(as.character(toPlot$variable))

    # REF052 vs Other
    p1 <- ggplot(toPlot, aes(x = variable, y = value, color = REF052)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
      geom_line(alpha = 0.65) +
      scale_color_manual("REF052", values = REF052_pal) +
      #coord_cartesian(ylim = c(0.75, 1.1)) +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm')
      ) +
      labs(
        y = "GC-corrected Fragment Midpoint Coverage",
        x = "Position relative to site (bp)",
        title = analysis
      )

    # Subtype
    p2 <- ggplot(toPlot, aes(x = variable, y = value, color = subtype)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
      geom_line(alpha = 0.65) +
      scale_color_manual("Subtype", values = cfDNA_subtype_pal) +
      #coord_cartesian(ylim = c(0.75, 1.1)) +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm')
      ) +
      labs(
        y = "GC-corrected Fragment Midpoint Coverage",
        x = "Position relative to site (bp)",
        title = analysis
      )

    filename <- paste0("data/results/figures/Misc/TF_Scoring/", analysis, "_coverage.png")
    png(filename, width=9, height=4, units='in', res = 600, pointsize=80)
    print(ggarrange(p1, p2, nrow = 1, ncol = 2))
    dev.off()
}

#-- modify from 5-cfDNA/ARCHE_CICADA.R/coverage_panel()

# get all griffin files
files <- list.files(
    paste0("data/rawdata/cfDNA/", dir),
    recursive = TRUE,
    pattern = "GC_corrected.coverage.tsv",
    full.names = TRUE
)
samples <- gsub("/.*", "", files)

# get covergage across positions
cov <- data.frame(matrix(nrow=0, ncol=134))
for (file in files) {
    df <- fread(file, data.table = FALSE)
    df <- df[, c("sample", "site_name", names(df)[1:132])] |> as.data.frame()
    cov <- rbind(cov, df)
}

# plot coverage plots
if (analysis == "UGT1A6") {
  # promoter
  toPlot <- cov[cov$site_name == "UGT1A6_promoter",]
  plot_coverage(toPlot, meta, "UGT1A6_promoter")
  # cCREs
  toPlot <- cov[cov$site_name == "UGT1A6_cCREs",]
  plot_coverage(toPlot, meta, "UGT1A6_cCREs")
} else {
  plot_coverage(cov, meta, analysis)
}
