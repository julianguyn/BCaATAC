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



###########################################################
# Compare TCGA scores with NMF
###########################################################

tumour_meta <- meta[meta$type == "tumour",]

scores <- zsc

scores <- as.data.frame(t(scores[,colnames(scores) %in% tumour_meta$sampleid]))
scores$ARCHE <- tumour_meta$arche[match(rownames(scores), tumour_meta$sampleid)]
scores$sampleID <- rownames(scores)
scores$rank <- as.character(tcga_meta$rank[match(scores$sampleID, gsub("\\.", "-", tcga_meta$Sample.Name))])


# --- heatmap

toPlot <- reshape2::melt(scores)
toPlot$variable <- factor(toPlot$variable, levels = paste0("ARCHE", 6:1))
toPlot <- toPlot[order(as.numeric(toPlot$rank)),]
toPlot$sampleID <- factor(toPlot$sampleID, levels = unique(toPlot$sampleID))

p1 <- ggplot(toPlot, aes(x = sampleID, y = "Assigned\nARCHE", fill = ARCHE)) +
    geom_tile() +
    scale_fill_manual(values = ARCHE_pal) +
    theme_void() +
    theme(legend.position = "none")

p2 <- ggplot(toPlot, aes(x = sampleID, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradientn("ARCHE\nExpression\nScore", colours = brewer.pal(9, "Blues")) +
    theme_void() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.text.y = element_text(size = 9, hjust=0)
    )

p <- p1 / p2 + plot_layout(heights = c(1, 7))

filename <- paste0("data/results/figures/Misc/all_scoring/", dir, "/tumour_heatmap.png")
ggsave(filename, p, w = 9, h = 4)

# --- individual waterfall plots

plot_ARCHE_waterfall <- function(arche, scores) { 

    toPlot <- scores[order(scores[[arche]], decreasing = TRUE),]
    toPlot$sampleID <- factor(toPlot$sampleID, levels = toPlot$sampleID)

    # get mismatch samples
    to_label <- c()
    n <- nrow(scores[scores$ARCHE == arche,])
    top_n <- toPlot[1:n,]
    if (nrow(top_n[top_n$ARCHE != arche,]) > 0) {
        ott <- toPlot[n+1:nrow(toPlot),]
        to_label <- c(top_n$sampleID[top_n$ARCHE != arche],  ott$sampleID[ott$ARCHE == arche])

        toPlot$label <- ifelse(toPlot$sampleID %in% to_label, "*", "")
    } else {
        toPlot$label <- ""
    }

    p <- ggplot(toPlot, aes(x = sampleID, y = .data[[arche]], fill = ARCHE)) +
        geom_bar(stat = "identity", color = "black") +
        geom_text(aes(label = label), vjust = -0.1, size = 6) +
        scale_fill_manual(values = ARCHE_pal) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
            axis.title.x = element_blank()
        )
    
    return(p)
}

compile_ARCHE_waterfall <- function(dir, scores) {

    p1 <- plot_ARCHE_waterfall("ARCHE1", scores)
    p2 <- plot_ARCHE_waterfall("ARCHE2", scores)
    p3 <- plot_ARCHE_waterfall("ARCHE3", scores)
    p4 <- plot_ARCHE_waterfall("ARCHE4", scores)
    p5 <- plot_ARCHE_waterfall("ARCHE5", scores)
    p6 <- plot_ARCHE_waterfall("ARCHE6", scores)

    p <- (p1 + p2) / (p3 + p4) / (p5 + p6) + plot_layout(guides="collect")

    filename <- paste0("data/results/figures/Misc/all_scoring/", dir, "/tumour_waterfall.png")
    ggsave(filename, p, w = 14, h = 10)

}

compile_ARCHE_waterfall(dir, scores)