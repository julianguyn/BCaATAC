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

missing <- c("REFLECT-0016-03", "REFLECT-0037-03", "REFLECT-0055-01", "REFLECT-0056-03")

###########################################################
# Plot heatmap
###########################################################

# helper function
plot_heatmap <- function(scores, group) {

  scores <- scores[scores$Subset == group,]

  df <- scores %>%
    select(Sample, Label, Score) %>%
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
      SeqBatch = meta$batch_num[match(colnames(df), meta$id_6b)],
      PipelineBatch = meta$Pipeline_batch[match(colnames(df), meta$id_6b)],
      col = list(Subtype = cfDNA_subtype_pal))

    filename <- paste0("data/results/figures/5-cfDNA/REFLECT/scores/", group, "_heatmap.png")
    png(filename, width = 9, height = 4, res = 600, units = "in")
    print(
        Heatmap(df, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
            row_names_gp = gpar(fontsize = 10), top_annotation = ha)
    )
    dev.off()
}

plot_heatmap(reflect, "10k")
plot_heatmap(reflect, "20k")
plot_heatmap(reflect, "50k")

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

plot_coverage <- function(cov, meta, group, subset = FALSE) {
    toPlot <- cov[grep(group, cov$site_name), ]
    toPlot <- reshape2::melt(toPlot)
    toPlot$subtype <- meta$Subtype_final[match(toPlot$sample, meta$id_6b)]
    toPlot$variable <- as.numeric(as.character(toPlot$variable))

    if (subset == TRUE) {
      toPlot$site_name <- sub("_.*", "", toPlot$site_name)
      toPlot <- toPlot[toPlot$subtype != "HER2",]
      toPlot$subtype[toPlot$subtype == "ER"] <- "ER+"
    }

    p <- ggplot(toPlot, aes(x = variable, y = value, color = subtype)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
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

    filename <- ifelse(
      subset == TRUE,
      paste0("data/results/figures/5-cfDNA/REFLECT/coverage/", group,"_subset_coverage.png"),
      paste0("data/results/figures/5-cfDNA/REFLECT/coverage/", group,"_coverage.png")
    )
    png(filename, width=8, height=5, units='in', res = 600, pointsize=80)
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

#-- plot coverage for figure
plot_coverage(cov, meta, "10k", subset = TRUE)
plot_coverage(cov, meta, "20k", subset = TRUE)
plot_coverage(cov, meta, "50k", subset = TRUE)

###########################################################
# Plot ranking by subtype
###########################################################

# helper function to plot tiles
plot_tiles <- function(df, arche) {

  toPlot <- df[df$ARCHE == arche,]
  toPlot <- toPlot[order(toPlot$Score, decreasing = TRUE),]
  toPlot$Rank <- 1:nrow(toPlot)

  p <- ggplot(toPlot, aes(x = Rank, y = 1, fill = Subtype)) +
    geom_tile() +
    scale_fill_manual(values = cfDNA_subtype_pal) +
    theme_void() +
    theme(
      axis.title.y = element_text(hjust = 1)
    ) +
    labs(y = arche)
  return(p)
}

# helper function to plot ranking tiles
plot_rank_scores <- function(group) {
  
  toPlot <- reflect[reflect$Subset == group,]
  
  a1 <- plot_tiles(toPlot, "ARCHE1")
  a2 <- plot_tiles(toPlot, "ARCHE2")
  a3 <- plot_tiles(toPlot, "ARCHE3")
  a4 <- plot_tiles(toPlot, "ARCHE4")
  a5 <- plot_tiles(toPlot, "ARCHE5")
  a6 <- plot_tiles(toPlot, "ARCHE6")

  p <- ggarrange(
    a1, a2, a3, a4, a5, a6,
    ncol = 1, nrow = 6,
    common.legend = TRUE
  )
  filename <- paste0("data/results/figures/5-cfDNA/REFLECT/rank/", group,"_Subtype.png")
  png(filename, width=5, height=3, units='in', res = 600, pointsize=80)
  print(p)
  dev.off()
}

plot_rank_scores("10k")
plot_rank_scores("20k")
plot_rank_scores("50k")

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
