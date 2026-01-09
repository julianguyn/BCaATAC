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
# Split into ARCHE subsets
###########################################################

c_unfiltered_10k <- c_unfiltered[c_unfiltered$Subset == "10k",]
c_unfiltered_20k <- c_unfiltered[c_unfiltered$Subset == "20k",]
c_unfiltered_50k <- c_unfiltered[c_unfiltered$Subset == "50k",]

c_bloodvser_10k <- c_bloodvser[c_bloodvser$Subset == "10k",]
c_bloodvser_20k <- c_bloodvser[c_bloodvser$Subset == "20k",]
c_bloodvser_50k <- c_bloodvser[c_bloodvser$Subset == "50k",]

###########################################################
# Plot correlation of ARCHE scores and TF
###########################################################

plot_ARCHE_TF_CICADA(c_unfiltered_10k, "CICADA-unfiltered_10k")
plot_ARCHE_TF_CICADA(c_unfiltered_20k, "CICADA-unfiltered_20k")
plot_ARCHE_TF_CICADA(c_unfiltered_50k, "CICADA-unfiltered_50k")

plot_ARCHE_TF_CICADA(c_bloodvser_10k, "CICADA-BloodvsER-25_10k")
plot_ARCHE_TF_CICADA(c_bloodvser_20k, "CICADA-BloodvsER-25_20k")
plot_ARCHE_TF_CICADA(c_bloodvser_50k, "CICADA-BloodvsER-25_50k")

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

c_unfiltered_10k <- make_toPlot(c_unfiltered_10k)
c_unfiltered_20k <- make_toPlot(c_unfiltered_20k)
c_unfiltered_50k <- make_toPlot(c_unfiltered_50k)

c_bloodvser_10k <- make_toPlot(c_bloodvser_10k)
c_bloodvser_20k <- make_toPlot(c_bloodvser_20k)
c_bloodvser_50k <- make_toPlot(c_bloodvser_50k)

###########################################################
# Plot ARCHE scores by TF group
###########################################################

plot_ARCHE_score_CICADA(c_unfiltered_10k, "CICADA-unfiltered_10k")
plot_ARCHE_score_CICADA(c_unfiltered_20k, "CICADA-unfiltered_20k")
plot_ARCHE_score_CICADA(c_unfiltered_50k, "CICADA-unfiltered_50k")

plot_ARCHE_score_CICADA(c_bloodvser_10k, "CICADA-BloodvsER-25_10k")
plot_ARCHE_score_CICADA(c_bloodvser_20k, "CICADA-BloodvsER-25_20k")
plot_ARCHE_score_CICADA(c_bloodvser_50k, "CICADA-BloodvsER-25_50k")

###########################################################
# Plot coverage plots
###########################################################

# helper function to plot coverage plots
coverage_panel <- function(folder, meta) {

  # get all griffin files
  dir <- paste0("data/rawdata/cfDNA/", folder)
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
  plot_coverage(cov, meta, "10k", folder)
  plot_coverage(cov, meta, "20k", folder)
  plot_coverage(cov, meta, "50k", folder)
}

coverage_panel("CICADA-unfiltered", meta)
coverage_panel("CICADA-BloodvsER-25", meta)