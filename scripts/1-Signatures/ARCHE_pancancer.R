# test ARCHE scoring across other cancer types

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

###########################################################
# Load in data
###########################################################

scores <- fread("data/rawdata/all_scoring/tcga_pancancer.Zscore.txt", data.table = FALSE)
rownames(scores) <- paste0("ARCHE", 1:6)
scores$V1 <- NULL

# load in pan-cancer metadata
meta <- read.csv("metadata/tcga_pancancer_meta.csv")

# load in BRCA metadata
tcga_meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
tcga_meta$ATAC.Seq.File.Name <- gsub("\\.", "-", tcga_meta$ATAC.Seq.File.Name)

###########################################################
# Map BRCA subtypes
###########################################################

meta$bca_subtype <- tcga_meta$Subtype[match(meta$sampleid, tcga_meta$ATAC.Seq.File.Name)]

###########################################################
# Plot heatmap
###########################################################

plot_heatmap <- function(scores, meta, label) {

    ha <- HeatmapAnnotation(
        Cancer = meta$cancer[match(colnames(scores), meta$sampleid)],
        Subtype = meta$bca_subtype[match(colnames(scores), meta$sampleid)],
        col = list(Cancer = cancer_type_pal, Subtype = subtype_pal),
        na_col = "white"
    )

    score_pal <- colorRamp2(seq(min(scores), max(scores), length = 3), c("#AD6A6C", "white", "#077293"))

    filename <- paste0("data/results/figures/1-Signatures/pancancer/", label, "_heatmap.png")
    cat("-----Saving plot to", filename, "\n")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(
        Heatmap(scores, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 3),
            row_names_gp = gpar(fontsize = 10), top_annotation = ha)
    )
    dev.off()

}

plot_heatmap(scores, meta, "zscores")