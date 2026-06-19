# test ARCHE scores with TCGA peaks as the bakground - plot heatmaps

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(matrixStats)
    library(ggh4x)
    library(reshape2)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(dplyr)
    library(readxl)
    library(data.table)
    library(patchwork)
    library(ComplexHeatmap)
    library(SummarizedExperiment)
    library(circlize)
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/plots/drug_response_ccls.R")
source("utils/palettes.R")
source("utils/bca_drugs.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/plots/drug_response_pdx.R")

# ---------------------------------------------------------
# Helper functions

get_arche_scores <- function(sample, filename, meta) {

    # load in arche scores
    scores <- read.table(filename)
    rownames(scores) <- paste0("ARCHE", 1:6)

    # standardize sample names of cell lines
    colnames(scores) <- sub("^X", "", gsub("\\.(?!$)", "-", colnames(scores), perl = TRUE))
    scores <- scores[, colnames(scores) %in% meta$filename]
    colnames(scores) <- meta$sampleid[match(colnames(scores), meta$filename)]
    scores <- scores[, order(colnames(scores))]
    if ("104987" %in% colnames(scores)) scores <- scores[, colnames(scores) != "104987"]
    if (sample == "pdxs") colnames(scores) <- map_pdx(colnames(scores))

    return(scores)
}

get_arche_devs <- function(sample, filename, meta) {

    scores <- fread(filename, data.table = FALSE) |> suppressWarnings()
    scores$V1 <- NULL
    rownames(scores) <- paste0("ARCHE", 1:6)

    # standardize sample names of cell lines
    colnames(scores) <- sub("^X", "", gsub("\\.(?!$)", "-", colnames(scores), perl = TRUE))
    scores <- scores[, colnames(scores) %in% meta$filename]
    colnames(scores) <- meta$sampleid[match(colnames(scores), meta$filename)]
    scores <- scores[, order(colnames(scores))]
    if ("104987" %in% colnames(scores)) scores <- scores[, colnames(scores) != "104987"]
    if (sample == "PDXs") colnames(scores) <- map_pdx(colnames(scores))

    return(scores)
}

###########################################################
# Prepare metadata
###########################################################

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")
dups <- meta$sampleid[duplicated(meta$sampleid)] # label dups
meta$sampleid[meta$sampleid %in% dups] <- paste0(meta$sampleid[meta$sampleid %in% dups], " (", meta$tech[meta$sampleid %in% dups], ")")

c_meta <- meta[meta$type == "cell_line", ]
p_meta <- meta[meta$type == "PDX", ]

###########################################################
# Load in cell line data
###########################################################

# load in arche zscores
zscore_cells <- get_arche_scores("cells", "data/rawdata/TCGA_bg/cells_50k.Zscore.txt", c_meta)

# load in arche deviations
deviat_cells <- get_arche_devs("cells", "data/rawdata/TCGA_bg/cells_50k.Deviations.txt", c_meta)

###########################################################
# Load in PDX data
###########################################################

# load in arche zscores
zscore_pdx <- get_arche_scores("pdxs", "data/rawdata/TCGA_bg/PDXs_50k.Zscore.txt", p_meta)

# load in arche deviations
deviat_pdx <- get_arche_devs("PDXs", "data/rawdata/TCGA_bg/PDXs_50k.Deviations.txt", p_meta)

###########################################################
# Try removing low deviation samples
###########################################################

# helper function to get sum of magnitude deviations
get_devs <- function(df, label, lim1, lim2) {

    devs <- colSums(abs(df)) |> as.data.frame()
    colnames(devs) <- "Sum"
    n1 <- paste0(as.character(length(devs[devs$Sum > lim1,])), "/", nrow(devs))
    n2 <- paste0(as.character(length(devs[devs$Sum > lim2,])), "/", nrow(devs))

    p <- ggplot(devs, aes(x = Sum)) +
        geom_histogram(fill = random_lightblue, color = "black", linewidth = 0.3) +
        geom_vline(xintercept = c(lim1, lim2), linetype = "dashed", color = "gray") +
        geom_text(stat = "bin",
            aes(label = after_stat(count)),
            vjust = -0.5, size = 3) +
        theme_bw() +
        ggtitle(paste0(label, ";  n >", lim1, ":", n1, ";  n >", lim2, ":", n2))
    filename <- paste0("data/results/figures/Misc/TCGA_bg/", label, ".png")
    ggsave(filename, p, w=5, h=4)

    to_keep <- rownames(devs[devs$Sum > lim1,,drop=FALSE])
    df <- df[,to_keep]
    return(df)
}

# get samples above threshold (cell lines)
zscore_cells_sumdev <- get_devs(zscore_cells, "zscore_cells", 55, 75)
deviat_cells_sumdev <- get_devs(deviat_cells, "deviat_cells", 0.25, 0.35)

# get samples above threshold (PDXs)
zscore_pdx_sumdev <- get_devs(zscore_pdx, "zscore_pdx", 55, 75)
deviat_pdx_sumdev <- get_devs(deviat_pdx, "deviat_pdx", 0.25, 0.4)

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

    filename <- paste0("data/results/figures/Misc/TCGA_bg_counts/", label, "_ARCHE_scores.png")
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

# plot heatmap of cell lines zscores (raw chromvar)
plot_ARCHE_scores_heatmap_counts(zscore_cells, "cells_zcores", c_meta)
plot_ARCHE_scores_heatmap_counts(zscore_cells, "cells_znorm_zscores", c_meta, znorm = TRUE)
plot_ARCHE_scores_heatmap_counts(zscore_cells_sumdev, "cells_zcores_sumdev", c_meta)
plot_ARCHE_scores_heatmap_counts(zscore_cells_sumdev, "cells_znorm_zscores_sumdev", c_meta, znorm = TRUE)

# heatmap of cell lines deviations
plot_ARCHE_scores_heatmap_counts(deviat_cells, "cells_deviations", c_meta)
plot_ARCHE_scores_heatmap_counts(deviat_cells, "cells_znorm_deviations", c_meta, znorm = TRUE)
plot_ARCHE_scores_heatmap_counts(deviat_cells_sumdev, "cells_deviations_sumdev", c_meta)
plot_ARCHE_scores_heatmap_counts(deviat_cells_sumdev, "cells_znorm_deviations_sumdev", c_meta, znorm = TRUE)


# plot heatmap of PDXs zscores (raw chromvar)
plot_ARCHE_scores_heatmap_counts(zscore_pdx, "pdxs_zcores", p_meta)
plot_ARCHE_scores_heatmap_counts(zscore_pdx, "pdxs_znorm_zscores", p_meta, znorm = TRUE)
plot_ARCHE_scores_heatmap_counts(zscore_pdx_sumdev, "pdxs_zcores_sumdev", p_meta)
plot_ARCHE_scores_heatmap_counts(zscore_pdx_sumdev, "pdxs_znorm_zscores_sumdev", p_meta, znorm = TRUE)

# heatmap of cell lines deviations
plot_ARCHE_scores_heatmap_counts(deviat_pdx, "pdxs_deviations", p_meta)
plot_ARCHE_scores_heatmap_counts(deviat_pdx, "pdxs_deviations_zscores", p_meta, znorm = TRUE)
plot_ARCHE_scores_heatmap_counts(deviat_pdx_sumdev, "pdxs_deviations_sumdev", p_meta)
plot_ARCHE_scores_heatmap_counts(deviat_pdx_sumdev, "pdxs_deviations_zscores_sumdev", p_meta, znorm = TRUE)
