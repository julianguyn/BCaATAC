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
cicada <- cicada[cicada$TF_group == "TF>10",]
cicada <- cicada[complete.cases(cicada),]

# get ARCHE scores for preclinical samples
preclinical <- score_arche_cfDNA("preclinical-ARCHE")
preclinical$ARCHE <- sub("sig", "ARCHE", preclinical$ARCHE)
preclinical$Sample <- sub("_merged", "", sub("_30", "", preclinical$Sample))
preclinical$Sample[preclinical$Sample == "CAMA1_mouse_ctDNA"] <- "CAMA1_xeno"

###########################################################
# Plot CICADA scores
###########################################################

# add labels
cicada$Sample_Type <- "Patient (CICADA)"
cicada$ARCHE <- factor(cicada$ARCHE, levels = paste0("ARCHE", 6:1))
cicada$plot_score <- log10((cicada$Score+0.01)*2) + 0.9
cicada$pheno <- cicada_meta$pheno_id[match(cicada$Sample, cicada_meta$sample_id)]
cicada$time <- cicada_meta$time_id[match(cicada$Sample, cicada_meta$sample_id)]

filename <- "data/results/figures/Misc/cicada_scores_time.png"
png(filename, width=6, height=4, units='in', res = 600, pointsize=80)
ggplot(cicada, aes(x = ARCHE, y = plot_score, fill = time)) +
    geom_boxplot() + 
    geom_jitter(aes(fill = time), shape = 21, stroke = 0.2,
                position = position_jitterdodge(jitter.width = 0.2)) +
    scale_fill_manual("Sample Type", values = c("#73937E", "#CEB992")) +
    scale_y_continuous(limits = c(-0.15, 0.85)) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        panel.border = element_rect(),
        legend.key.size = unit(0.7, 'cm')
    ) +
    labs(y = "ARCHE Score")
dev.off()

###########################################################
# Plot preclinical scores
###########################################################

preclinical$Sample <- factor(
    preclinical$Sample,
    levels = c("CAMA1_xeno", "CAMA1", "MCF7", "BPTO95"),
    labels = c("Xenograft (CAMA-1)", "CCL (CAMA-1)", "CCL (MCF-7)", "Organoid"))

filename <- "data/results/figures/Misc/preclinical_scores_time.png"
png(filename, width=6, height=4, units='in', res = 600, pointsize=80)
ggplot(preclinical, aes(x = ARCHE, y = Score, fill = Sample)) + 
    geom_bar(stat = "identity", position = position_dodge(), color = "black") + 
    geom_hline(yintercept = 0) + 
    scale_fill_manual(
        "Sample Type",
        values = c("#655560", "#C57B57", "#F1AB86", "#F7DBA7"),
        labels = ) +
    scale_y_continuous(limits = c(-0.03, 0.55)) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        panel.border = element_rect(),
        legend.key.size = unit(0.7, 'cm')
    ) + 
    labs(y = "ARCHE Score")
dev.off()