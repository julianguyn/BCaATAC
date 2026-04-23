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

# ---------------------------------------------------------
# sumdev
zscore_sumdev <- switch(
    samples,
    cells_50k = 150,
    PDXs_50k = 100
)
deviat_sumdev <- switch(
    samples,
    cells_50k = 0.5,
    PDXs_50k = 0.3
)

# helper function to subset for low deviation samples
subset_sumdev <- function(df, lim) {
    devs <- colSums(abs(df)) |> as.data.frame()
    colnames(devs) <- "Sum"
    to_keep <- rownames(devs[devs$Sum > lim,,drop=FALSE])
    cat("Keeping", length(to_keep), "samples\n")
    df <- df[,to_keep]
}

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
if ("104987" %in% colnames(arche_dev)) arche_dev <- arche_dev[,-which(colnames(arche_dev) == "104987")]

###########################################################
# Plot heatmaps
###########################################################

# helper function to make heatmap
plot_ARCHE_scores_heatmap_counts <- function(df, label, meta, znorm = FALSE, subset_dev = FALSE, lim = NULL) {

    # normalize
    if (znorm == TRUE) {
        cat("Normalizing\n")
        toPlot <- znorm(df)
        toPlot <- toPlot[, colSums(is.na(toPlot)) == 0]
        df <- df[,colnames(df) %in% colnames(toPlot)]
    } else {
        toPlot <- df
    }
    
    score_pal <- colorRamp2(seq(min(toPlot), max(toPlot), length = 3), c("#AD6A6C", "white", "#077293"))
    count_pal <- colorRamp2(seq(min(colSums(abs(df))), max(colSums(abs(df))), length = 3), c("#DFDFDF", "#989898", "#202020"))

    ha <- HeatmapAnnotation(
        Subtype = meta[match(colnames(df), meta$sampleid),]$subtype,
        Tech = meta[match(colnames(df), meta$sampleid),]$tech,
        SumDevs = colSums(abs(df)),
        col = list(Subtype = subtype_pal, Tech = tech_pal, SumDevs = count_pal))

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
plot_ARCHE_scores_heatmap_counts(subset_sumdev(arche_zscores, zscore_sumdev), paste0(samples, "_zscores_sumdev"), meta)
plot_ARCHE_scores_heatmap_counts(subset_sumdev(arche_zscores, zscore_sumdev), paste0(samples, "_znorm_zscores_sumdev") , meta, znorm = TRUE)

# plot heatmap of deviation scores
plot_ARCHE_scores_heatmap_counts(arche_dev, paste0(samples, "_deviations"), meta)
plot_ARCHE_scores_heatmap_counts(arche_dev, paste0(samples, "_znorm_deviations") , meta, znorm = TRUE)
plot_ARCHE_scores_heatmap_counts(subset_sumdev(arche_dev, deviat_sumdev), paste0(samples, "_deviations_sumdev"), meta)
plot_ARCHE_scores_heatmap_counts(subset_sumdev(arche_dev, deviat_sumdev), paste0(samples, "_znorm_deviations_sumdev") , meta, znorm = TRUE)

###########################################################
# Compare rank of ARCHE scoring
###########################################################

order_samples <- colnames(arche_zscores) # already ordered

# helper function to prepare data for merge
format_arche_scores <- function(df, label, norm = FALSE) {
    to_keep <- order_samples[order_samples %in% colnames(df)]
    df <- df[,match(to_keep, colnames(df))]
    if (norm == TRUE) df <- znorm(df)
    df$ARCHE <- paste0("ARCHE", 1:6)
    to_merge <- reshape2::melt(df)
    to_merge$Label <- label
    return(to_merge)
}


toPlot <- rbind(
    format_arche_scores(arche_zscores, "ZScores"),
    format_arche_scores(arche_dev, "Deviations"),
    format_arche_scores(arche_zscores, "znormZScores", norm = TRUE),
    format_arche_scores(arche_dev, "znormDeviations", norm = TRUE),
    format_arche_scores(subset_sumdev(arche_zscores, zscore_sumdev), "Zscores_sumdev"),
    format_arche_scores(subset_sumdev(arche_dev, deviat_sumdev),"Deviations_sumdev" ),
    format_arche_scores(subset_sumdev(arche_zscores, zscore_sumdev), "znormZScores_sumdev", norm = TRUE),
    format_arche_scores(subset_sumdev(arche_dev, deviat_sumdev), "znormDeviations_sumdev", norm = TRUE)
)

plot_order <- function(arche, label) {
    subset <- toPlot[toPlot$Label == label & toPlot$ARCHE == arche,]
    cell_order <- as.character(rev(subset$variable[order(subset$value)]))
    subset$variable <- factor(subset$variable, levels = cell_order)

    subtype <- meta[match(subset$variable[order(subset$value)], meta$sampleid),]
    subtype$sampleid <- factor(subtype$sampleid, levels = cell_order)

    p <- ggplot(subtype, aes(x = sampleid, y = arche, fill = subtype)) +
        geom_tile(color = "black") +
        scale_fill_manual("Subtype", values = subtype_pal) +
        theme_void() +
        theme(
            axis.text.y = element_text(hjust = 0.5),
            legend.position = "none"
        ) 

    return(p)
}

scores_list <- c(
    "ZScores", "Deviations", "znormZScores", "znormDeviations",
    "ZScores_sumdev", "Deviations_sumdev", "znormZScores_sumdev", "znormDeviations_sumdev"
)
for (score_type in scores_list) {
    p1 <- plot_order("ARCHE1", score_type)
    p2 <- plot_order("ARCHE2", score_type)
    p3 <- plot_order("ARCHE3", score_type)
    p4 <- plot_order("ARCHE4", score_type)
    p5 <- plot_order("ARCHE5", score_type)
    p6 <- plot_order("ARCHE6", score_type)

    p <- (p1 / p2 / p3 / p4 / p5 / p6) +
    plot_annotation(title = score_type, theme = theme(plot.title = element_text(hjust = 0.5)))

    filename <- paste0("data/results/figures/Misc/ARCHE_counts/", samples, "_", score_type, "_rank.png")
    ggsave(filename, p, width = 7, height = 2)
}


###########################################################
# Old code: Comparison znorm scores across ARCHE scoring
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
    
    if (sub("_.*", "", samples) == "cells") {
        ggsave(filename, p1+p2+plot_layout(heights = c(1, 2)), width = 9, height = 2.25)
    } else {
        ggsave(filename, p1+p2+plot_layout(heights = c(1, 2)), width = 10, height = 2)
    }
}

#for (arche in paste0("ARCHE", 1:6)) {
#    plot_ordered_arche_scores(arche, "ZScores")
#    plot_ordered_arche_scores(arche, "Deviations")
#}