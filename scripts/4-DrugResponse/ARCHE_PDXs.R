# load libraries
suppressPackageStartupMessages({
    library(readxl)
    library(ggplot2)
    library(RColorBrewer)
    library(dplyr)
    library(ROCR)
    library(ggpubr)
    library(ggh4x)
    library(grid)
    library(gridExtra)
    library(matrixStats)
})

source("utils/plots/drug_response_pdx.R")
source("utils/plots/drug_response_pdx_indivplots.R")
source("utils/palettes.R")
source("utils/get_data.R")
source("utils/mappings.R")
source("utils/plots/ARCHE_scores_heatmap.R")

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

p_meta <- meta[meta$type == "PDX", ]

# load in arche scores
pdxs_20k <- get_arche_scores("pdxs", "k20", p_meta)
pdxs_50k <- get_arche_scores("pdxs", "k50", p_meta)
pdxs_all <- get_arche_scores("pdxs", "all", p_meta)

# pdx drug response data
xeva <- get_xeva("full")

###########################################################
# Normalize ARCHE scores
###########################################################

norm_20k <- znorm(pdxs_20k)
norm_50k <- znorm(pdxs_50k)
norm_all <- znorm(pdxs_all)

###########################################################
# Assign ARCHE scores
###########################################################

# helper function
format_ARCHE <- function(xeva, scores) {
    scores <- t(scores) |> as.data.frame()
    df <- xeva[xeva$patient.id %in% rownames(scores),]
    df$ARCHE6 <- df$ARCHE5 <- df$ARCHE4 <- df$ARCHE3 <- df$ARCHE2 <- df$ARCHE1 <- NA
    for (i in 1:nrow(df)) {
        sample <- df$patient.id[i]
        df$ARCHE1[i] <- scores[rownames(scores) == sample,]$ARCHE1
        df$ARCHE2[i] <- scores[rownames(scores) == sample,]$ARCHE2
        df$ARCHE3[i] <- scores[rownames(scores) == sample,]$ARCHE3
        df$ARCHE4[i] <- scores[rownames(scores) == sample,]$ARCHE4
        df$ARCHE5[i] <- scores[rownames(scores) == sample,]$ARCHE5
        df$ARCHE6[i] <- scores[rownames(scores) == sample,]$ARCHE6
        
    }
    return(df)
}

# get ARCHE scores
xeva_20k <- format_ARCHE(xeva, pdxs_20k)
xeva_50k <- format_ARCHE(xeva, pdxs_50k)
xeva_all <- format_ARCHE(xeva, pdxs_all)

xeva_norm_20k <- format_ARCHE(xeva, norm_20k)
xeva_norm_50k <- format_ARCHE(xeva, norm_50k)
xeva_norm_all <- format_ARCHE(xeva, norm_all)

###########################################################
# Assess ARCHE drug response associations
###########################################################

dir <- "rmNergizKomal-fullXEVA"

# create the data dir
path <- paste0("data/results/data/4-DrugResponse/PDX/", dir)
if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
}
# create the figure dirs
to_create <- c("PDX20k", "PDX20K_norm", "PDX50k", "PDX50k_norm", "PDXall", "PDXall_norm")
for (subfolder in to_create) {
    path <- paste0("data/results/figures/4-DrugResponse/PDX/", dir, "/", subfolder)
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
}


x1 <- assess_ARCHE_PDX(xeva_20k, "PDX20k", dir) |> suppressWarnings()
x2 <- assess_ARCHE_PDX(xeva_50k, "PDX50k", dir) |> suppressWarnings()
x3 <- assess_ARCHE_PDX(xeva_all, "PDXall", dir) |> suppressWarnings()

n1 <- assess_ARCHE_PDX(xeva_norm_20k, "PDX20k_norm", dir) |> suppressWarnings()
n2 <- assess_ARCHE_PDX(xeva_norm_50k, "PDX50k_norm", dir) |> suppressWarnings()
n3 <- assess_ARCHE_PDX(xeva_norm_all, "PDXall_norm", dir) |> suppressWarnings()

##todo: why does it keep printing NULL out

###########################################################
# Subset for significant associations
###########################################################

# compile results
compile <- rbind(x1, x2, x3, n1, n2, n3)

# set thresholds
pc <- 0.4
pval <- 0.1

# plot bubble plots
plot_PDX_bubbles(compile, "BR", dir)
plot_PDX_bubbles(compile, "BAR", dir)