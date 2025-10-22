# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(reshape2)
    library(ComplexHeatmap)
    library(circlize)
})

source("source/DataExploration/helper.R")
source("source/DataExploration/plots.R")
source("source/DrugResponse/helper.R")
source("source/palettes.R")

###########################################################
# Load in data
###########################################################

# get gene counts from PSets
ubr2 <- load_bca_RNA()
gray <- get_pset_rna("DrugResponse/data/PSet_GRAY2017.rds")
gcsi <- get_pset_rna("DrugResponse/data/gCSI.rds")
ccle <- get_pset_rna("DrugResponse/data/CCLE.rds")

###########################################################
# Melt for correlations
###########################################################

# common genes
genes <- intersect(rownames(ubr2), rownames(ccle))

# melt
ubr2_m <- reshape2::melt(ubr2[match(genes, rownames(ubr2)),order(colnames(ubr2))])
gray_m <- reshape2::melt(gray[match(genes, rownames(gray)),order(colnames(gray))])
gcsi_m <- reshape2::melt(gcsi[match(genes, rownames(gcsi)),order(colnames(gcsi))])
ccle_m <- reshape2::melt(ccle[match(genes, rownames(ccle)),order(colnames(ccle))])


###########################################################
# Correlate RNA-seq expression
###########################################################

p1 <- corr_pset_rna(ubr2_m, gray_m)
p2 <- corr_pset_rna(ubr2_m, gcsi_m)
p3 <- corr_pset_rna(ubr2_m, ccle_m)
p4 <- corr_pset_rna(gray_m, gcsi_m)
p5 <- corr_pset_rna(gray_m, ccle_m)
p6 <- corr_pset_rna(gcsi_m, ccle_m)

###########################################################
# Plot RNA-seq correlation matrix
###########################################################

plot_rna_corr(p1, p2, p3, p4, p5, p6)