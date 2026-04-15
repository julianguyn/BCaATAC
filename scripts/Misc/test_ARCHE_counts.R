# test ARCHE counts as signature score

suppressPackageStartupMessages({
    library(data.table)
    library(matrixStats)
    library(ComplexHeatmap)
    library(SummarizedExperiment)
    library(circlize)
    library(reshape2)
    library(ggplot2)
    library(patchwork)
})

set.seed(101)

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/palettes.R")
source("utils/plots/ARCHE_scores_heatmap.R")

# set samples to load in
args <- commandArgs(trailingOnly = TRUE)
samples <- args[1]
valid <- c("cells_20k", "cells_50k", "cells_all", "PDXs_20k", "PDXs_50k", "PDXs_all")
if (is.na(samples) || !samples %in% valid) {
  stop(
    sprintf("Invalid analysis argument '%s'. Must be one of: %s",
            samples, paste(valid, collapse = ", ")),
    call. = FALSE
  )
}

arche_thres <- switch(
    sub(".*_", "group_", samples),
    group_20k = "k20",
    group_50k = "k50",
    group_all = "all"
)

###########################################################
# Load in data
###########################################################

cat("Loading in data\n")

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# label dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta$sampleid[meta$sampleid %in% dups] <- paste0(meta$sampleid[meta$sampleid %in% dups], " (", meta$tech[meta$sampleid %in% dups], ")")

# load in zscore deviation scoring (original)
arche_zscores <- get_arche_scores(tolower(sub("_.*", "", samples)), arche_thres, meta)

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
    #cat("Number of sites:", length(sites), "\n")
    subset <- counts[sites,]
    #print(dim(subset))
    scores <- colSums(subset)
    #print(scores)
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
    cat("-----Saving plot to", filename, "\n")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(
        Heatmap(toPlot, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
            row_names_gp = gpar(fontsize = 10), top_annotation = ha)
    )
    dev.off()
}

cat("Making plots\n")

# plot heatmap of zscore deviations (chromvar)
plot_ARCHE_scores_heatmap_counts(arche_zscores, paste0(samples, "_zscores"), meta)
plot_ARCHE_scores_heatmap_counts(arche_zscores, paste0(samples, "_znorm_zscores") , meta, znorm = TRUE)

# plot heatmap of deviation scores
plot_ARCHE_scores_heatmap_counts(arche_dev, paste0(samples, "_deviations"), meta)
plot_ARCHE_scores_heatmap_counts(arche_dev, paste0(samples, "_znorm_deviations") , meta, znorm = TRUE)

# plot heatmap of count matrix scores
plot_ARCHE_scores_heatmap_counts(arche_scores, paste0(samples, "_countmat"), meta)
plot_ARCHE_scores_heatmap_counts(arche_scores, paste0(samples, "_znorm_countmat") , meta, znorm = TRUE)

###########################################################
# Comparison znorm scores across ARCHE scoring
###########################################################

order_samples <- colnames(arche_zscores) # already ordered

# helper function to prepare data for merge
format_arche_scores <- function(df, label) {
    df <- df[,match(order_samples, colnames(df))]
    df <- znorm(df)
    df$ARCHE <- paste0("ARCHE", 1:6)
    to_merge <- reshape2::melt(df)
    to_merge$Label <- label
    return(to_merge)
}

toPlot <- rbind(
    format_arche_scores(arche_zscores, "ZScores"),
    format_arche_scores(arche_dev, "Deviations")
    #format_arche_scores(arche_scores, "Counts") - not representative
)

plot_ordered_arche_scores <- function(arche, order) {

    subset <- toPlot[toPlot$ARCHE == arche,]
    to_order <- subset[subset$Label == order,]
    cell_order <- rev(to_order$variable[order(to_order$value)])

    subset$variable <- factor(subset$variable, levels = cell_order)
    if (order == "Deviations") subset$Label <- factor(subset$Label, levels = c("ZScores", "Deviations"))

    subtype <- meta[match(to_order$variable[order(to_order$value)], meta$sampleid),]
    subtype$sampleid <- factor(subtype$sampleid, levels = cell_order)

    p1 <- ggplot(subtype, aes(x = sampleid, y = "Subtype ", fill = subtype)) +
        geom_tile(color = "black") +
        scale_fill_manual("Subtype", values = subtype_pal) +
        theme_void() +
        theme(
            axis.text.y = element_text(),
            legend.position = "none"
        ) 
    
    p2 <- ggplot(subset, aes(x = variable, y = Label, fill = value)) +
        geom_tile(color = "black") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.key.size = unit(0.3, 'cm')
        ) +
        scale_fill_gradient2(
            paste0(arche, "\nScore"),
            low = "#B05C58", mid = "#F8F1F8", high = "#077293", limits = c(-2, 2)
        )

    filename <- paste0("data/results/figures/Misc/ARCHE_counts/", samples, "_", arche, "_", order, ".png")
    ggsave(filename, p1+p2+plot_layout(heights = c(1, 2)), width = 9, height = 2.25)
}

for (arche in paste0("ARCHE", 1:6)) {
    plot_ordered_arche_scores(arche, "ZScores")
    plot_ordered_arche_scores(arche, "Deviations")
}
