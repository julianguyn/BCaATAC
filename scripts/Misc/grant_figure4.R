# temp script, ideally move these into correct scripts afterwards

# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(viridis)
  library(data.table)
  library(dplyr)
  library(circlize)
})

source("utils/score_arche_cfDNA.R")
source("utils/plots/cfdna.R")
source("utils/palettes.R")
source("utils/mappings.R")

###########################################################
# Load in data
###########################################################

# load in CICADA metadata
cicada_meta <- read.csv("data/rawdata/cfDNA/CICADA-ARCHE/sample_metadata.csv")

# get ARCHE scores for cicada
cicada <- score_arche_cfDNA("CICADA-BloodvsER-25", meta = cicada_meta)
arche_keep <- c(paste0("ARCHE", c(1, 3:6), "_50k"), "ARCHE2_20k")
cicada <- cicada[cicada$Label %in% arche_keep,]
cicada$TF_group <- ifelse(cicada$TF > 5, "TF>5", "Not")
cicada <- cicada[cicada$TF_group == "TF>5",]

# get ARCHE scores for preclinical samples
preclinical <- score_arche_cfDNA("preclinical-ARCHE")
preclinical$ARCHE <- sub("sig", "ARCHE", preclinical$ARCHE)
preclinical$Sample <- sub("_merged", "", sub("_30", "", preclinical$Sample))
preclinical$Sample[preclinical$Sample == "CAMA1_mouse_ctDNA"] <- "CAMA1_xeno"

###########################################################
# Plot preclinical scores
###########################################################

preclinical$Sample <- factor(
    preclinical$Sample,
    levels = c("CAMA1_xeno", "CAMA1", "MCF7", "BPTO95"),
    labels = c("Xenograft (CAMA-1)", "CCL (CAMA-1)", "CCL (MCF-7)", "Organoid"))

filename <- "data/results/figures/Misc/preclinical_scores_time.png"

p1 <- ggplot(preclinical, aes(x = ARCHE, y = Score, fill = Sample)) + 
    geom_bar(stat = "identity", position = position_dodge(), color = "black") + 
    geom_hline(yintercept = 0) + 
    scale_fill_manual(
        "Sample Type",
        values = c("#655560", "#C57B57", "#F1AB86", "#F7DBA7"),
        labels = ) +
    scale_y_continuous(limits = c(-0.01, 0.55), expand = c(0,0)) +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(color = "black", linewidth = 0.2), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    ) + 
    labs(y = "ARCHE Score")
p <- p1 / p2 + plot_layout(height = c(20, 1))
ggsave(filename, p, width = 4, height = 3.5)

###########################################################
# Plot CICADA scores
###########################################################

# add labels
cicada$Sample_Type <- "Patient (CICADA)"
cicada$ARCHE <- factor(cicada$ARCHE, levels = paste0("ARCHE", 1:6))
cicada$pheno <- cicada_meta$pheno_id[match(cicada$Sample, cicada_meta$sample_id)]
cicada$time <- cicada_meta$time_id[match(cicada$Sample, cicada_meta$sample_id)]
cicada <- cicada[complete.cases(cicada),]

filename <- "data/results/figures/Misc/cicada_scores_time.png"

p1 <- ggplot(cicada, aes(x = ARCHE, y = Score, fill = time)) +
    geom_boxplot() + 
    geom_jitter(aes(fill = time), shape = 21, stroke = 0.2,
                position = position_jitterdodge(jitter.width = 0.2)) +
    scale_fill_manual("Sample Type", values = c("#73937E", "#CEB992")) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.background = element_rect(color = "black", linewidth = 0.2), 
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    ) +
    labs(y = "ARCHE Score")

p2 <- ggplot(cicada, aes(x = ARCHE, y = "", fill = ARCHE)) +  
    geom_tile(color = "black") +
    scale_fill_manual(values = ARCHE_pal) +
    theme_void() +
    theme(
      axis.text.x = element_text(size = 9, vjust = 0),
      legend.position = "none",
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    )
p <- p1 / p2 + plot_layout(height = c(20, 1))
ggsave(filename, p, width = 4, height = 3.5)