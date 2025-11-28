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
# Load in data
###########################################################

# load in metadata
dir <- "data/rawdata/cfDNA/CICADA-ARCHE"
meta <- read.csv(paste0(dir, "/sample_metadata.csv"))

###########################################################
# Get ARCHE scores
###########################################################

c_unfiltered <- score_arche_cfDNA("CICADA-unfiltered", meta)
c_bloodvser <- score_arche_cfDNA("CICADA-BloodvsER-25", meta)

###########################################################
# Plot correlation of ARCHE scores and TF
###########################################################

plot_ARCHE_TF_CICADA(c_unfiltered, "CICADA-unfiltered")
plot_ARCHE_TF_CICADA(c_bloodvser, "CICADA-BloodvsER-25")

###########################################################
# Get ARCHE score summary stats
###########################################################

# helper function to compile results
make_toPlot <- function(scores) {

  toPlot <- rbind(
    summary_ARCHE(scores, "All Samples"),
    summary_ARCHE(scores[scores$TF_group == "TF>10",], "TF>10"),
    summary_ARCHE(scores[scores$TF_group == "TF<10",], "TF<10")
  )

  # formating df
  toPlot$TF_group <- factor(toPlot$TF_group, levels = c("All Samples", "TF>10", "TF<10"))
  toPlot$label <- factor(toPlot$label, levels = c("Baseline", "Progression", "On treatment control"))
  toPlot <- toPlot[complete.cases(toPlot),]

  return(toPlot)
}

c_unfiltered <- make_toPlot(c_unfiltered)
c_bloodvser <- make_toPlot(c_bloodvser)

###########################################################
# Plot ARCHE scores by TF group
###########################################################

plot_ARCHE_score_CICADA(c_unfiltered, "CICADA-unfiltered")
plot_ARCHE_score_CICADA(c_bloodvser, "CICADA-BloodvsER-25")