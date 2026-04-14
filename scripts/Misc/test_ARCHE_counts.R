# test ARCHE counts as signature score

suppressPackageStartupMessages({
    library(data.table)
    library(matrixStats)
    library(ComplexHeatmap)
    library(SummarizedExperiment)
    library(circlize)
    library(reshape2)
    library(ggplot2)
})

set.seed(101)

source("utils/get_data.R")
source("utils/palettes.R")
source("utils/plots/ARCHE_scores_heatmap.R")

# set samples to load in
samples <- "cells_20k"

###########################################################
# Load in data
###########################################################

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# label dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta$sampleid[meta$sampleid %in% dups] <- paste0(meta$sampleid[meta$sampleid %in% dups], " (", meta$tech[meta$sampleid %in% dups], ")")

# load in sample ARCHE deviation scores
arche_dev <- read.table(paste0("data/rawdata/ARCHE_counts/", samples, ".Deviations.txt"))
colnames(arche_dev) <- sub("^X", "", gsub("\\.(?!$)", "-", colnames(arche_dev), perl = TRUE))
colnames(arche_dev) <- meta$sampleid[match(colnames(arche_dev), meta$filename)]

# load in ARCHE counts
counts <- paste0("data/rawdata/ARCHE_counts/", samples, ".counts_filtered.Rdata")
load(counts)
counts <- assay(counts_filtered, "counts")
counts <- as.data.frame.matrix(counts)

# load in ARCHE consensus matrix
arche_mat <- readRDS(paste0("data/rawdata/ARCHE_counts/", samples, ".consensus.Binarymat.repeats.rds"))
rownames(arche_mat) <- paste(arche_mat$seqnames, arche_mat$start, arche_mat$end, sep = "_")
arche_mat <- arche_mat[,-which(colnames(arche_mat) %in% c("seqnames", "start", "end"))]

###########################################################
# Get ARCHE scores from counts and consensus matrices
###########################################################

arche_scores <- data.frame(matrix(nrow=0, ncol=ncol(counts)))
for (arche in colnames(arche_mat)) {
    cat("--------Starting", arche, "\n")
    sites <- rownames(arche_mat)[arche_mat[[arche]] == 1]
    cat("Number of sites:", length(sites), "\n")
    subset <- counts[sites,]
    print(dim(subset))
    scores <- colSums(subset)
    print(scores)
    arche_scores <- rbind(arche_scores, scores)
}
rownames(arche_scores) <- colnames(arche_mat)
colnames(arche_scores) <- colnames(counts)
colnames(arche_scores) <- meta$sampleid[match(colnames(arche_scores), meta$filename)]

###########################################################
# Plot heatmaps
###########################################################

# helper function to make heatmap
plot_ARCHE_scores_heatmap_counts <- function(df, label, meta, znorm = FALSE) {

    # normalize
    if (znorm == TRUE) {
        cat("Normalizing\n")
        toPlot <- znorm(df)
    } else {
        toPlot <- df
    }
    
    score_pal <- colorRamp2(seq(min(toPlot), max(toPlot), length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))
    count_pal <- colorRamp2(seq(min(colSums(df)), max(colSums(df)), length = 3), c("#E8D6CB", "#AD6A6C", "#5D2E46"))

    ha <- HeatmapAnnotation(
        Subtype = meta[match(colnames(df), meta$sampleid),]$subtype,
        Tech = meta[match(colnames(df), meta$sampleid),]$tech,
        PeakCount = colSums(df),
        col = list(Subtype = subtype_pal, Tech = tech_pal, PeakCount = count_pal))

    filename <- paste0("data/results/figures/Misc/ARCHE_counts/", label, "_ARCHE_scores.png")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(
        Heatmap(toPlot, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
            row_names_gp = gpar(fontsize = 10), top_annotation = ha)
    )
    dev.off()
}

# plot heatmap of deviation scores
plot_ARCHE_scores_heatmap_counts(arche_dev, paste0(samples, "_deviations"), meta)
plot_ARCHE_scores_heatmap_counts(arche_dev, paste0(samples, "_znorm_deviations") , meta, znorm = TRUE)

# plot heatmap of count matrix scores
plot_ARCHE_scores_heatmap_counts(arche_scores, paste0(samples, "_countmat"), meta)
plot_ARCHE_scores_heatmap_counts(arche_scores, paste0(samples, "_znorm_countmat") , meta, znorm = TRUE)
