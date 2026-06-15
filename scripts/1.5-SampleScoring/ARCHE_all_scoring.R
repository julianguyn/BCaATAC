# test ARCHE scores with all samples + TCGA (cell + pdx + tcga)

# load libraries
suppressPackageStartupMessages({
    library(matrixStats)
    library(data.table)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    library(patchwork)
})

set.seed(101)

source("utils/palettes.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/plots/test_ARCHE_all_scoring.R")

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]

valid <- c("cell_pdx_tcga", "cell_pdx_pdo_tcga")
if (is.na(analysis) || !analysis %in% valid) {
  stop(
    sprintf("Invalid analysis argument '%s'. Must be one of: %s",
            analysis, paste(valid, collapse = ", ")),
    call. = FALSE
  )
}

# ---------------------------------------------------------
# Helper functions

get_scores <- function(filename, meta) {
    scores <- fread(filename, data.table = FALSE) |> suppressWarnings()
    scores$V1 <- NULL
    rownames(scores) <- paste0("ARCHE", 1:6)
    colnames(scores) <- meta$sampleid[match(sub("_peaks.*", "", colnames(scores)), meta$filename)]
    if ("104987" %in% colnames(scores)) scores <- scores[, colnames(scores) != "104987"]
    return(scores)
}

###########################################################
# Prepare metadata
###########################################################

# read in sample metadata
sample_meta <- read.csv("metadata/lupien_metadata.csv")
dups <- sample_meta$sampleid[duplicated(sample_meta$sampleid)] # label dups
sample_meta$sampleid[sample_meta$sampleid %in% dups] <- paste0(sample_meta$sampleid[sample_meta$sampleid %in% dups], " (", sample_meta$tech[sample_meta$sampleid %in% dups], ")")

# load in TCGA metadata
tcga_meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
tcga_meta$Sample.Name[tcga_meta$Sample.Name == "TCGA.A2.A0T4"] <- paste0("TCGA.A2.A0T4-", 1:2)

# load in PDO metadata
pdo_meta <- read.csv("metadata/pdo_metadata.csv")

# compile metadata
if (dir == "cell_pdx_tcga") {
    meta <- data.frame(
        filename = c(sample_meta$filename, gsub("\\.", "-", tcga_meta$ATAC.Seq.File.Name)),
        sampleid = c(sample_meta$sampleid, gsub("\\.", "-", tcga_meta$Sample.Name)),
        type = c(sample_meta$type, rep("tumour", nrow(tcga_meta))),
        subtype = c(sample_meta$subtype, tcga_meta$Subtype),
        arche = c(rep("NA", nrow(sample_meta)), tcga_meta$ARCHE)
    )
} else if (dir == "cell_pdx_pdo_tcga") {
    meta <- data.frame(
        filename = c(sample_meta$filename, pdo_meta$filename, gsub("\\.", "-", tcga_meta$ATAC.Seq.File.Name)),
        sampleid = c(sample_meta$sampleid, pdo_meta$sampleid, gsub("\\.", "-", tcga_meta$Sample.Name)),
        type = c(sample_meta$type, pdo_meta$type, rep("tumour", nrow(tcga_meta))),
        subtype = c(sample_meta$subtype, pdo_meta$subtype, tcga_meta$Subtype),
        arche = c(rep("NA", nrow(sample_meta)+nrow(pdo_meta)), tcga_meta$ARCHE)
    )
}

# condense subtypes
meta$subtype[meta$subtype == "unknown"] <- "Not Available"
meta$subtype[meta$subtype == "Basal"] <- "Basal/TNBC"
meta$subtype[meta$subtype == "TNBC"] <- "Basal/TNBC"
meta$subtype[meta$subtype == "LumA"] <- "Luminal/ER+"
meta$subtype[meta$subtype == "LumB"] <- "Luminal/ER+"
meta$subtype[meta$subtype == "ER"] <- "Luminal/ER+"
meta$subtype[meta$subtype == "Her2"] <- "HER2"

###########################################################
# Load in scores
###########################################################

# zscores
zsc <- get_scores(paste0("data/rawdata/all_scoring/", dir, ".Zscore.txt"), meta)

# load in arche deviations
dev <- get_scores(paste0("data/rawdata/all_scoring/", dir, ".Deviations.txt"), meta)

###########################################################
# Remove low deviation samples
###########################################################

# get samples above threshold
zscore_sumdev <- get_devs(zsc, dir, meta, "zscore", 100)
deviat_sumdev <- get_devs(dev, dir, meta, "deviations", 0.25)

###########################################################
# Plot heatmaps
###########################################################

# plot heatmap of zscores
plot_ARCHE_scores_heatmap_counts(zsc, dir, "zcores", meta)
plot_ARCHE_scores_heatmap_counts(zsc, dir, "znorm_zscores", meta, znorm = TRUE)

# heatmap of deviations
plot_ARCHE_scores_heatmap_counts(dev, dir, "deviations", meta)
plot_ARCHE_scores_heatmap_counts(dev, dir, "znorm_deviations", meta, znorm = TRUE)

###########################################################
# PCA to check alignment
###########################################################

# pca outfiles generated from 1.5-SampleScoring/ARCHE_all_PCA.R
load(paste0("data/rawdata/all_scoring/", dir, ".RData"))

plot_pca(pca_res, dir, meta)

###########################################################
# Compare rank of ARCHE scoring
###########################################################

for (arche in paste0("ARCHE", 1:6)) {

    p1 <- plot_ARCHE_rank(zsc, arche, "ZScores")
    p2 <- plot_ARCHE_rank(zscore_sumdev, arche, "ZScores-sumdev")
    p3 <- plot_ARCHE_rank(zsc, arche, "ZScores (znorm)", norm = TRUE)
    p4 <- plot_ARCHE_rank(zscore_sumdev, arche, "ZScores-sumdev (znorm)", norm = TRUE)

    p5 <- plot_ARCHE_rank(dev, arche, "Deviations")
    p6 <- plot_ARCHE_rank(deviat_sumdev, arche, "Deviations-sumdev")
    p7 <- plot_ARCHE_rank(dev, arche, "Deviations (znorm)", norm = TRUE)
    p8 <- plot_ARCHE_rank(deviat_sumdev, arche, "Deviations-sumdev (znorm)", norm = TRUE)


    p <- p1 / p2 / p3 / p4 / p5 / p6 / p7 / p8
    filename <- paste0("data/results/figures/Misc/all_scoring/", dir, "/rank_scores_", arche, ".png")
    ggsave(filename, p, width = 13, height = 6)
}


###########################################################
# Compare TCGA scores with NMF
###########################################################

tumour_meta <- meta[meta$type == "tumour",]

plot_tumour_ARCHE_scores(zsc, dir, tumour_meta, "ZScores")
plot_tumour_ARCHE_scores(znorm(zsc), dir, tumour_meta, "ZScores (znorm)")
plot_tumour_ARCHE_scores(dev, dir, tumour_meta, "Deviations")
plot_tumour_ARCHE_scores(znorm(dev), dir, tumour_meta, "Deviations (znorm)")