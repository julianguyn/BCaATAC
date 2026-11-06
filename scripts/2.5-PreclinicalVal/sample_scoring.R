# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(matrixStats)
    library(reshape2)
    library(patchwork)
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/palettes.R")
source("utils/plots/ARCHE_scores_heatmap.R")

# ---------------------------------------------------------
# Cell lines

###########################################################
# Prepare metadata
###########################################################

# read in sample metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]

# get cell lines
c_meta <- meta[meta$type == "cell_line", ]

###########################################################
# Load in cell line data
###########################################################

# zscores
zscore_cells <- get_arche_scores("data/rawdata/all_scoring/cell_tcga.Zscore.txt", c_meta)
normzs_cells <- znorm(zscore_cells)

# sumdevs
zscore_cells_sumdev <- get_arche_sumdevs(zscore_cells, "zscore_cells", plot = TRUE)
normzs_cells_sumdev <- znorm(zscore_cells_sumdev)

###########################################################
# Plot subtype enrichment scores
###########################################################

plot_subtype_enrichment(zscore_cells, c_meta, "zscore_cells")
plot_subtype_enrichment(normzs_cells, c_meta, "normzs_cells")
plot_subtype_enrichment(zscore_cells_sumdev, c_meta, "zscore_cells_sumdev")
plot_subtype_enrichment(normzs_cells_sumdev, c_meta, "normzs_cells_sumdev")


# ---------------------------------------------------------
# PDXs

###########################################################
# Prepare metadata
###########################################################

# read in sample metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]

# remove second dups
dups <- meta$sampleid[duplicated(meta$sampleid)] # only PDXs have dups
meta_t <- meta[!(meta$sampleid %in% dups & meta$tech == "tina"), ]
meta_k <- meta[!(meta$sampleid %in% dups & meta$tech == "komal"), ]

###########################################################
# Load in PDX data
###########################################################

p_meta_t <- meta_t[meta_t$type == "PDX", ]
p_meta_k <- meta_k[meta_k$type == "PDX", ]

# load in arche zscores and normalize
zscore_t <- get_arche_scores("data/rawdata/all_scoring/pdx_tcga.Zscore.txt", p_meta_t)
normzs_t <- znorm(zscore_t)
zscore_k <- get_arche_scores("data/rawdata/all_scoring/pdx_tcga.Zscore.txt", p_meta_k)
normzs_k <- znorm(zscore_k)

###########################################################
# Plot subtype enrichment scores
###########################################################

plot_subtype_enrichment(zscore_t, meta_t, "pdx_zscore_cells")
plot_subtype_enrichment(normzs_t, meta_t, "pdx_normzs_cells")
plot_subtype_enrichment(zscore_k, meta_k, "pdx_zscore_cells_sumdev")
plot_subtype_enrichment(normzs_k, meta_k, "pdx_normzs_cells_sumdev")