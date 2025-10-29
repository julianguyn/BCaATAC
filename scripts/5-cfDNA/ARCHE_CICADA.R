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

dir <- "data/rawdata/cfDNA/CICADA-ARCHE"

# load in metadata
meta <- read.csv(paste0(dir, "/sample_metadata.csv"))

# samples to remove from analysis (did not finish running)
to_rm <- c("CICADA-0031-01", "CICADA-0011-02")
meta <- meta[-which(meta$sample_id %in% to_rm), ]

###########################################################
# Get ARCHE scores
###########################################################

scores <- score_arche_cfDNA(dir)

###########################################################
# Get variables from metadata
###########################################################

scores$label <- meta$time_id[match(scores$Sample, meta$sample_id)]
scores$TF <- meta$metrics_tf[match(scores$Sample, meta$sample_id)]
scores$pheno <- meta$pheno_id[match(scores$Sample, meta$sample_id)]

scores$TF_group <- ifelse(scores$TF > 10, "TF>10", "TF<10")

###########################################################
# Get ARCHE score summary stats
###########################################################

# helper function to compute mean and SE of ARCHE scores

summary_ARCHE <- function(df, TF_group) {
    summary_df <- df %>%
        group_by(ARCHE, label) %>%
        summarise(
            mean_score = mean(Score, na.rm = TRUE),
            se_score = sd(Score, na.rm = TRUE) / sqrt(n())
        ) %>%
        ungroup()
    summary_df <- as.data.frame(summary_df)
    summary_df$TF_group <- TF_group
    return(summary_df)
}

toPlot <- rbind(
  summary_ARCHE(scores, "All Samples"),
  summary_ARCHE(scores[scores$TF_group == "TF>10",], "TF>10"),
  summary_ARCHE(scores[scores$TF_group == "TF<10",], "TF<10")
)

###########################################################
# Plot ARCHE scores
###########################################################

# formating df
toPlot$TF_group <- factor(toPlot$TF_group, levels = c("All Samples", "TF>10", "TF<10"))
toPlot$label <- factor(toPlot$label, levels = c("Baseline", "Progression", "On treatment control"))
toPlot <- toPlot[complete.cases(toPlot),]

# plot ARCHE scores by TF group
plot_ARCHE_score_CICADA(toPlot)

# plot correlation between TF and ARCHE scores
plot_ARCHE_TF_CICADA(scores)
