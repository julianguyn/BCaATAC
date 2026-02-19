# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(readxl)
  library(ggplot2)
  library(data.table)
  library(tidyverse)
  library(umap)
})

source("utils/score_arche_cfDNA.R")
source("utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

dir <- "data/rawdata/cfDNA/REFLECT-unfiltered"

# load in metadata
meta <- read_excel(paste0(dir, "/REFLECT-6B_metadata.xlsx"), sheet = 1) |> as.data.frame()

###########################################################
# Format metadata
###########################################################

meta$Subtype_final[meta$Subtype_final == "ER+/HER2-"] <- "ER+"
meta$Subtype_final[meta$Subtype_final == "TBNC"] <- "TNBC"
meta$id_6b <- gsub("_", "-", gsub("_LB.*", "", meta$library_id))

meta$Subtype_final[meta$Subtype_final == "ER+"] <- "ER"
meta$Subtype_final[meta$Subtype_final == "HER2+"] <- "HER2"

###########################################################
# Get ARCHE scores
###########################################################

reflect <- score_arche_cfDNA("REFLECT-unfiltered")
reflect$Subtype <- meta$Subtype_final[match(reflect$Sample, meta$id_6b)]
reflect$SeqBatch <- meta$batch_num[match(reflect$Sample, meta$id_6b)]
reflect$PipelineBatch <- meta$Pipeline_batch[match(reflect$Sample, meta$id_6b)]

# check incomplete files
all_samples <- list.files(dir)
table(all_samples %in% reflect$Sample)
all_samples[-which(all_samples %in% reflect$Sample)]

missing <- c(
  "REFLECT-0005-03","REFLECT-0007-03","REFLECT-0016-03",
  "REFLECT-0021-03", "REFLECT-0025-03", "REFLECT-0027-01",
  "REFLECT-0037-03", "REFLECT-0043-01", "REFLECT-0044-01",
  "REFLECT-0050-03", "REFLECT-0055-01", "REFLECT-0056-03",
  "REFLECT-0062-03")

###########################################################
# Plot scores
###########################################################

# helper function
plot_scores <- function(scores, group, metric = "cc") {

  toPlot <- scores[scores$Subset == group,]

  # get metric
  label <- switch(
    metric,
    cc = "Central Covarege",
    mc = "Mean Coverage"
  )

  p <- ggplot(toPlot, aes(x = ARCHE, y = Score, fill = Subtype)) +
    geom_boxplot() +
    geom_jitter(
        shape = 21,
        position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
    ) +
    scale_fill_manual(values = cfDNA_subtype_pal) +
    ylim(-0.05, 0.65) +
    theme_minimal() +
    theme(
      panel.border = element_rect(),
      legend.key.size = unit(0.5, 'cm')
    ) + labs(y = paste0("Score (1 - ", label, ")"))

  filename <- paste0("data/results/figures/5-cfDNA/REFLECT/scores/", group, "_", metric, "_scores.png")
  png(filename, width=6, height=4, units='in', res = 600, pointsize=80)
  print(p)
  dev.off()
}

plot_scores(reflect, "10k")
plot_scores(reflect, "20k")
plot_scores(reflect, "50k")

###########################################################
# Plot coverage plots
# modify from utils/plots/cfdna.R
###########################################################

plot_coverage <- function(cov, meta, group) {
    toPlot <- cov[grep(group, cov$site_name), ]
    toPlot <- reshape2::melt(toPlot)
    toPlot$subtype <- meta$Subtype_final[match(toPlot$sample, meta$id_6b)]
    toPlot$variable <- as.numeric(as.character(toPlot$variable))

    p <- ggplot(toPlot, aes(x = variable, y = value, color = subtype)) +
      geom_line(alpha = 0.65) +
      facet_wrap(~ site_name, scales = "fixed", nrow = 2) +
      scale_color_manual("Subtype", values = cfDNA_subtype_pal) +
      coord_cartesian(ylim = c(0.3, 1.1)) +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm')
      ) +
      labs(
        y = "GC-corrected Fragment Midpoint Coverage",
        x = "Position relative to site (bp)"
      )

    png(paste0("data/results/figures/5-cfDNA/REFLECT/coverage/", group,"_coverage.png"), width=9, height=6, units='in', res = 600, pointsize=80)
    print(p)
    dev.off()
}

#-- modify from 5-cfDNA/ARCHE_CICADA.R/coverage_panel()

# get all griffin files
files <- list.files(
    dir,
    recursive = TRUE,
    pattern = "GC_corrected.coverage.tsv"
)
samples <- gsub("/.*", "", files)

# get covergage across positions
cov <- data.frame(matrix(nrow=0, ncol=134))
for (file in files) {
    df <- fread(paste0(dir, "/", file), data.table = FALSE)
    df <- df[, c("sample", "site_name", names(df)[1:132])] |> as.data.frame()
    cov <- rbind(cov, df)
}

# plot coverage plots
plot_coverage(cov, meta, "10k")
plot_coverage(cov, meta, "20k")
plot_coverage(cov, meta, "50k")

###########################################################
# Plot batch effects
###########################################################

plot_batch_effects <- function(group) {
  toPlot <- reflect %>%
    filter(Subset == group) %>%
    select(Sample, ARCHE, Score) %>%
    pivot_wider(
      names_from = ARCHE,
      values_from = Score
    ) %>%
    column_to_rownames(var = "Sample")

  # perform umap
  umap_scores <- umap(toPlot)$layout
  colnames(umap_scores) <- paste0("UMAP", 1:2)

  # get annotations
  anno <- unique(reflect[reflect$Subset == group,c("Sample", "Subtype", "SeqBatch", "PipelineBatch")])
  anno <- anno[match(rownames(umap_scores), anno$Sample),]
  umap_scores <- cbind(umap_scores, anno) |> as.data.frame()
  umap_scores$PipelineBatch <- factor(umap_scores$PipelineBatch)
  umap_scores$SeqBatch <- factor(umap_scores$SeqBatch)

  filename <- paste0("data/results/figures/5-cfDNA/REFLECT/batches/", group,"_Subtype.png")
  png(filename, width=4, height=3, units='in', res = 600, pointsize=80)
  print(ggplot(umap_scores, aes(x = UMAP1, y = UMAP2, fill = Subtype)) +
    geom_point(shape = 21, size = 2) +
    scale_fill_manual(values = cfDNA_subtype_pal) +
    theme_minimal() +
    theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.key.size = unit(0.5, 'cm')
        ))
  dev.off()

  filename <- paste0("data/results/figures/5-cfDNA/REFLECT/batches/", group,"_PipelineBatch.png")
  png(filename, width=4, height=3, units='in', res = 600, pointsize=80)
  print(ggplot(umap_scores, aes(x = UMAP1, y = UMAP2, fill = PipelineBatch)) +
    geom_point(shape = 21, size = 2) +
    #scale_fill_manual(values = cfDNA_subtype_pal) +
    theme_minimal() +
    theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.key.size = unit(0.5, 'cm')
        ))
  dev.off()

  filename <- paste0("data/results/figures/5-cfDNA/REFLECT/batches/", group,"_SeqBatch.png")
  png(filename, width=4, height=3, units='in', res = 600, pointsize=80)
  print(ggplot(umap_scores, aes(x = UMAP1, y = UMAP2, fill = SeqBatch)) +
    geom_point(shape = 21, size = 2) +
    #scale_fill_manual(values = cfDNA_subtype_pal) +
    theme_minimal() +
    theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          legend.key.size = unit(0.5, 'cm')
        ))
  dev.off()
}

plot_batch_effects("10k")
plot_batch_effects("20k")
plot_batch_effects("50k")
