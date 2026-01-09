suppressPackageStartupMessages({
    library(data.table)
    library(matrixStats)
    library(ComplexHeatmap)
    library(circlize)
    library(reshape2)
    library(ggplot2)
})

set.seed(101)

source("utils/get_data.R")
source("utils/palettes.R")
source("utils/plots/ARCHE_scores_heatmap.R")

###########################################################
# Load in data
###########################################################

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# label dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta$sampleid[meta$sampleid %in% dups] <- paste0(meta$sampleid[meta$sampleid %in% dups], " (", meta$tech[meta$sampleid %in% dups], ")")

c_meta <- meta[meta$type == "cell_line", ]
p_meta <- meta[meta$type == "PDX", ]

# load in arche scores
cells_20k <- get_arche_scores("cells", "k20", c_meta)
cells_50k <- get_arche_scores("cells", "k50", c_meta)
cells_all <- get_arche_scores("cells", "all", c_meta)

pdxs_20k <- get_arche_scores("pdxs", "k20", p_meta)
pdxs_50k <- get_arche_scores("pdxs", "k50", p_meta)
pdxs_all <- get_arche_scores("pdxs", "all", p_meta)

###########################################################
# Plot ARCHE scores
###########################################################

plot_ARCHE_scores_heatmap(cells_20k, "cells_20k", c_meta)
plot_ARCHE_scores_heatmap(cells_50k, "cells_50k", c_meta)
plot_ARCHE_scores_heatmap(cells_all, "cells_all", c_meta)

plot_ARCHE_scores_heatmap(pdxs_20k, "pdxs_20k", p_meta)
plot_ARCHE_scores_heatmap(pdxs_50k, "pdxs_50k", p_meta)
plot_ARCHE_scores_heatmap(pdxs_all, "pdxs_all", p_meta)

###########################################################
# Compare ARCHE signature subsets
###########################################################

plot_ARCHE_scores_compare(cells_20k, cells_50k, cells_all, "cells")
plot_ARCHE_scores_compare(pdxs_20k, pdxs_50k, pdxs_all, "pdxs")

# ---------------------------------------------------------
# Remove duplicates: nergiz and komal

###########################################################
# Load in data
###########################################################

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]

# remove komal dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "komal"), ]

c_meta <- meta[meta$type == "cell_line", ]
p_meta <- meta[meta$type == "PDX", ]

# load in arche scores
cells_20k <- get_arche_scores("cells", "k20", c_meta)
cells_50k <- get_arche_scores("cells", "k50", c_meta)
cells_all <- get_arche_scores("cells", "all", c_meta)

pdxs_20k <- get_arche_scores("pdxs", "k20", p_meta)
pdxs_50k <- get_arche_scores("pdxs", "k50", p_meta)
pdxs_all <- get_arche_scores("pdxs", "all", p_meta)

###########################################################
# Plot ARCHE scores
###########################################################

plot_ARCHE_scores_heatmap(cells_20k, "cells_20k", c_meta)
plot_ARCHE_scores_heatmap(cells_50k, "cells_50k", c_meta)
plot_ARCHE_scores_heatmap(cells_all, "cells_all", c_meta)

plot_ARCHE_scores_heatmap(pdxs_20k, "pdxs_20k", p_meta)
plot_ARCHE_scores_heatmap(pdxs_50k, "pdxs_50k", p_meta)
plot_ARCHE_scores_heatmap(pdxs_all, "pdxs_all", p_meta)

###########################################################
# Compare ARCHE signature subsets
###########################################################

plot_ARCHE_scores_compare(cells_20k, cells_50k, cells_all, "cells")
plot_ARCHE_scores_compare(pdxs_20k, pdxs_50k, pdxs_all, "pdxs")