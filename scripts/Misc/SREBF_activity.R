# SREBF1/2 binding site accessibility for Mitchell

suppressPackageStartupMessages({
    library(data.table)
    library(matrixStats)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
})

source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# label dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta$sampleid[meta$sampleid %in% dups] <- paste0(meta$sampleid[meta$sampleid %in% dups], " (", meta$tech[meta$sampleid %in% dups], ")")

meta <- meta[meta$type == "cell_line", ]

# chromvar scoring 
scores <- read.table("data/rawdata/tracks/cells_SREBF.Zscore.txt")

###########################################################
# Format scores (from get_arche_scores())
###########################################################

# standardize sample names of cell lines
colnames(scores) <- sub("^X", "", gsub("\\.(?!$)", "-", colnames(scores), perl = TRUE))
scores <- scores[, colnames(scores) %in% meta$filename]
colnames(scores) <- meta$sampleid[match(colnames(scores), meta$filename)]
scores <- scores[, order(colnames(scores))]

###########################################################
# Plot output of chromvar (from plots/ARCHE_scores_heatmap.R)
###########################################################

#' Modified ARCHE_scores_heatmap()
#' @param scores data.frame. Deviations from chromvar
#' @param label string. Label for figure filename
#' @param h int. Image height for figure
#' 
TF_scores_heatmap <- function(scores, label, h) {

    # set colours for plotting
    lim <- max(c(abs(min(scores)), max(scores)))
    score_pal <- colorRamp2(seq(-lim, lim, length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))

    ha <- HeatmapAnnotation(
        Subtype = meta[match(colnames(scores), meta$sampleid),]$subtype,
        col = list(Subtype = subtype_pal)
    )

    filename <- paste0("data/results/figures/Misc/", label, ".png")
    png(filename, width = 10, height = h, res = 600, units = "in")
    print(Heatmap(scores, cluster_rows = FALSE, name = paste0(label, "\nBinding\nScore"), col = score_pal,
        column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 10), top_annotation = ha))
    dev.off()
}

# plot across all sites
TF_scores_heatmap(scores, "SREBF", 3)

# plot for individual sites
for (i in 1:nrow(scores)) {
    TF_scores_heatmap(scores[i,], rownames(scores)[i], 2.75)
}
