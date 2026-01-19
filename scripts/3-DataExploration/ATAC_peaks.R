# load libraries
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
})

source("utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

# load in binary matrices of consensus peaks
df_cells <- readRDS("data/rawdata/ccls/cells_50k.consensus.Binarymat.rds")
df_pdxs <- readRDS("data/rawdata/pdx/PDXs_50k.consensus.Binarymat.rds")

# read in metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# label dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta$dup <- ifelse(meta$sampleid %in% dups, meta$sampleid, "unique")
meta$sampleid[meta$sampleid %in% dups] <- paste0(meta$sampleid[meta$sampleid %in% dups], " (", meta$tech[meta$sampleid %in% dups], ")")

###########################################################
# Format peaks for correlation
###########################################################

# helper function to format dataframes
format_mats <- function(df, meta) {

    # format peak names
    rownames(df) <- paste(df$seqnames, df$start, df$end, sep = ":")
    df$seqnames <- df$start <- df$end <- NULL

    # format sample names
    colnames(df) <- meta$sampleid[match(colnames(df), meta$filename)]
    return(df)
}

df_cells <- format_mats(df_cells, meta)
df_pdxs <- format_mats(df_pdxs, meta)


###########################################################
# Perform correlations
###########################################################

corr_cells <- cor(df_cells, use = "complete.obs")
corr_pdxs <- cor(df_pdxs, use = "complete.obs")

###########################################################
# Perform correlations
###########################################################

# helper function to plot heatmap
plot_heatmap <- function(corr_df, label) {

    # set colours for plotting
    score_pal <- colorRamp2(seq(min(corr_df), max(corr_df), length = 3), c("#BBCDE5", "#5189BD", "#49516F"))

    # colours for labeling dups
    m <- meta[meta$sampleid %in% colnames(corr_df),]
    l <- length(unique(m$dup)) - 1
    pal <- c(brewer.pal(brewer.pal.info["Paired", "maxcolors"], "Paired"), "#BAB9B8", "#6E6E6E")
    pal <- c(pal[1:l], "unique" = "white")
    names(pal)[1:l] <- unique(m$dup)[-which(unique(m$dup) == "unique")]

    # set subtype and clustering order annotation
    col_fun <- colorRamp2(c(0, 30, 50), c("#CEF0F2", "#42B6A3", "#446D73"))
    ha <- HeatmapAnnotation(
        Subtype = meta$subtype[match(colnames(corr_df), meta$sampleid)],
        Tech = meta$tech[match(colnames(corr_df), meta$sampleid)],
        Dup = meta$dup[match(colnames(corr_df), meta$sampleid)],
        col = list(Subtype = subtype_pal, Tech = tech_pal, Dup = pal)
    )

    # plot heatmap
    ht <- Heatmap(
        corr_df,
        name = "Correlation",
        col = score_pal,
        column_title = "Tumour Sample",
        column_title_side = "bottom",
        show_column_names = TRUE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        top_annotation = ha
    )
    ht <- draw(ht)

    filename <- paste0("data/results/figures/3-DataExploration/ATACheatmaps/", label, "_atac_heatmap.png")
    png(filename, width = 10, height = 9, res = 600, units = "in")
    print(draw(ht))
    dev.off()
}

plot_heatmap(corr_cells, "cells")
plot_heatmap(corr_pdxs, "pdxs")