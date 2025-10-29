# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(data.table)
  library(dplyr)
})

source("utils/score_arche_cfDNA.R")
source("utils/plots/cfdna.R")
source("utils/palettes.R")
source("utils/mappings.R")

###########################################################
# Get ARCHE scores
###########################################################

dir <- "data/rawdata/cfDNA/preclinical-ARCHE"
scores <- score_arche_cfDNA(dir)

###########################################################
# Map sample names
###########################################################

scores$Sample <- factor(
  scores$Sample, 
  levels = c("CAMA1_30", "MCF7_merged", "BPTO95", "CAMA1_mouse_ctDNA_merged"),
  labels = c("CCL (CAMA-1)", "CCL (MCF-7)", "Organoid", "Xenograft")
)

###########################################################
# Plot ARCHE scores
###########################################################

plot_ARCHE_score_preclinical(scores)
