setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(data.table)
  library(dplyr)
})

source("source/cfDNA/helper.R")
source("source/palettes.R")

###########################################################
# Specify files to load in
###########################################################

# get list of files to open
path <- "cfDNA/data/2025-09-09_julia_signature_results/"

# load in metadata
meta <- read.csv(paste0(path, "sample_metadata.csv"))

# samples to remove from analysis (did not finish running)
to_rm <- c("CICADA-0031-01", "CICADA-0011-02")
meta <- meta[-which(meta$sample_id %in% to_rm), ]

# get sample names for each group
ontrtctl <- meta$sample_id[meta$time_id == "On treatment control"]
baseline <- meta$sample_id[meta$time_id == "Baseline"]
progress <- meta$sample_id[meta$time_id == "Progression"]

# get files for each group
b_files <- paste0(path, baseline, "/", baseline, ".GC_corrected.coverage.tsv")
p_files <- paste0(path, progress, "/", progress, ".GC_corrected.coverage.tsv")
c_files <- paste0(path, ontrtctl, "/", ontrtctl, ".GC_corrected.coverage.tsv")

###########################################################
# Get ARCHE scores
###########################################################

baseline <- score_ARCHE(b_files, "Baseline")
progress <- score_ARCHE(p_files, "Progression")
ontrtctl <- score_ARCHE(c_files, "On treatment control")

###########################################################
# Get tumour fraction values
###########################################################

baseline$TF <- meta$metrics_tf[match(baseline$Sample, meta$sample_id)]
progress$TF <- meta$metrics_tf[match(progress$Sample, meta$sample_id)]
ontrtctl$TF <- meta$metrics_tf[match(ontrtctl$Sample, meta$sample_id)]

###########################################################
# Get ARCHE score summary stats
###########################################################

# get ARCHE score summary stats
df <- rbind(baseline, progress, ontrtctl)
df$TF_group <- ifelse(df$TF > 10, "TF>10", "TF<10")
toPlot <- rbind(
  summary_ARCHE(df, "All Samples"),
  summary_ARCHE(df[df$TF_group == "TF>10",], "TF>10"),
  summary_ARCHE(df[df$TF_group == "TF<10",], "TF<10")
)

# formating df
toPlot$TF_group <- factor(toPlot$TF_group, levels = c("All Samples", "TF>10", "TF<10"))
toPlot$label <- factor(toPlot$label, levels = c("Baseline", "Progression", "On treatment control"))
toPlot <- toPlot[complete.cases(toPlot),]

###########################################################
# Plot ARCHE scores by TF group
###########################################################

plot_ARCHE_score_CICADA(toPlot)