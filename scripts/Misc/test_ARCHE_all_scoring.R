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
})

set.seed(101)

source("utils/palettes.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/plots/test_ARCHE_all_scoring.R")

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

# compile metadata
meta <- data.frame(
    filename = c(sample_meta$filename, gsub("\\.", "-", tcga_meta$ATAC.Seq.File.Name)),
    sampleid = c(sample_meta$sampleid, gsub("\\.", "-", tcga_meta$Sample.Name)),
    type = c(sample_meta$type, rep("tumour", nrow(tcga_meta))),
    subtype = c(sample_meta$subtype, tcga_meta$Subtype),
    arche = c(rep("NA", nrow(sample_meta)), tcga_meta$ARCHE)
)

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
zsc <- get_scores("data/rawdata/all_scoring/cell_pdx_tcga.Zscore.txt", meta)

# load in arche deviations
dev <- get_scores("data/rawdata/all_scoring/cell_pdx_tcga.Deviations.txt", meta)

###########################################################
# Remove low deviation samples
###########################################################

# get samples above threshold
zscore_sumdev <- get_devs(zsc, meta, "zscore", 100)
deviat_sumdev <- get_devs(dev, meta, "deviations", 0.25)

###########################################################
# Plot heatmaps
###########################################################

# plot heatmap of zscores
plot_ARCHE_scores_heatmap_counts(zsc, "zcores", meta)
plot_ARCHE_scores_heatmap_counts(zsc, "znorm_zscores", meta, znorm = TRUE)

# heatmap of deviations
plot_ARCHE_scores_heatmap_counts(dev, "deviations", meta)
plot_ARCHE_scores_heatmap_counts(dev, "znorm_deviations", meta, znorm = TRUE)

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

ggplot(toPlot, aes(x = sampleID, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradientn("ARCHE\nExpression\nScore", colours = brewer.pal(9, "Blues")) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
        axis.title.x = element_blank()
    )

# --- individual waterfall plots

for (arche in paste0("ARCHE", 1:6)) {

    toPlot <- scores[order(scores[[arche]], decreasing = TRUE),]
    toPlot$sampleID <- factor(toPlot$sampleID, levels = toPlot$sampleID)

    ggplot(toPlot, aes(x = sampleID, y = .data[[arche]], fill = ARCHE)) +
        geom_bar(stat = "identity", color = "black") +
        scale_fill_manual(values = ARCHE_pal) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
            axis.title.x = element_blank()
        )

}