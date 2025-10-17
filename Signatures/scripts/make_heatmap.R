setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(data.table)
    library(ggplot2)
    library(reshape2)
    library(dplyr)
    library(RColorBrewer)
    library(ggh4x)
    library(ggpubr)
    library(grid)
    library(gridExtra)
})

source("source/Signatures/plots.R")
source("source/palettes.R")

###########################################################
# Load in data
###########################################################

# read in meta data
meta <- read.csv("Signatures/data/TCGA_subtype_label.csv")

# load in matrix file from NMF
mat <- read.table("Signatures/results/data/ATAC_heatmap_rank6.png.order.matrix", header = T)

###########################################################
# Format signature data for plotting
###########################################################

# format for plotting and get subtype
mat$Signature <- paste0("ARCHE", 1:6)
mat$Signature <- factor(mat$Signature, levels = paste0("ARCHE",6:1))
mat <- reshape2::melt(mat)
mat$subtype <- meta$bca_subtype[match(gsub("X", "", mat$variable), gsub("-", "\\.", meta$File.Name))]

# save signature assignment
mat$signature_assign <- ""
for (sample in mat$variable) {
    tmp <- mat[mat$variable == sample,]
    mat[mat$variable == sample,]$signature_assign <- as.character(tmp[which.max(tmp$value),]$Signature)
}

###########################################################
# Plot ATAC-Signature heatmap
###########################################################

plot_ARCHE_heatmap(mat)
