suppressPackageStartupMessages({
    library(data.table)
    library(matrixStats)
    library(ComplexHeatmap)
    library(circlize)
    library(reshape2)
})

source("utils/get_data.R")
source("utils/palettes.R")
sourcE("utils/plots/ARCHE_scores_heatmap.R")

###########################################################
# Load in data
###########################################################

# read in scoring on 20k sites
p2 <- read.table("data/rawdata/pdx/chromvar-Nov132025/BCa_PDXs.Zscore.txt")
colnames(p2) <- sub("_.*", "", sub("X", "", colnames(p2)))
p2 <- p2[,-which(colnames(p2) %in% "104987")]


# read in scoring on all sites (old)
pa <- read.table("data/rawdata/pdx/chromvar/bca_sign.Zscore.txt")
colnames(pa) <- sub("X", "", colnames(pa))
rownames(pa) <- paste0("ARCHE", 1:6)

# read in scoring on all sites (new)
p3 <- read.table("data/rawdata/pdx/chromvar-Nov132025-all/BCa_PDXs.Zscore.txt")
colnames(p3) <- sub("_.*", "", sub("X", "", colnames(p3)))
# note: 104987 has scores here

# read in scoring on cell lines (old)
c1 <- get_arche_cells()

# get cell line subtypes
c_subtypes <- get_cell_subtype()

###########################################################
# Plot ARCHE scores
###########################################################

plot_ARCHE_scores_heatmap(p2, "PDX_20k")
plot_ARCHE_scores_heatmap(pa, "PDX_all_old")
plot_ARCHE_scores_heatmap(p3, "PDX_all_new")
plot_ARCHE_scores_heatmap(c1, "CCLs_all_old", c_subtypes)

###########################################################
# Compare 20k to all ARCHE sites
###########################################################