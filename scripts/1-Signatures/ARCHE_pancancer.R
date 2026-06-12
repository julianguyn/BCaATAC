# test ARCHE scoring across other cancer types

# load libraries
suppressPackageStartupMessages({
    library(matrixStats)
    library(data.table)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
})



x <- fread("data/rawdata/all_scoring/tcga_pancancer.Zscore.txt")

# load in metadata
meta <- read.csv("metadata/tcga_pancancer_meta.csv")
