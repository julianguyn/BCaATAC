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
# Get cfDNA ARCHE scores
###########################################################

dir <- "preclinical-ARCHE-2026"
scores <- score_arche_cfDNA(dir)
scores$ARCHE <- sub("sig", "ARCHE", scores$ARCHE)
scores$Sample <- sub("_merged", "", sub("_30", "", scores$Sample))
scores$Sample[scores$Sample == "CAMA1_mouse_ctDNA"] <- "CAMA1_xeno"

###########################################################
# Plot ARCHE scores
###########################################################

plot_ARCHE_score_preclinical(scores, "2026")
plot_ARCHE_score_preclinical(scores, "2026")
plot_ARCHE_score_preclinical(scores, "2026")