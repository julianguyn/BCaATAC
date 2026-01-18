# format tracks from Chip-Atlas for chromVar + plots for Mitchell

suppressPackageStartupMessages({
    library(data.table)
    library(matrixStats)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
})

source("utils/palettes.R")

###########################################################
# Format TF binding sites bed file
###########################################################

# helper function to format bed files from Chip-Atlas
format_bed <- function(qval) {
    # load in and format file
    filename <- paste0("data/rawdata/tracks/Oth.Brs.", qval, ".NFE2L2.AllCell.bed")
    bed <- fread(filename, skip = "track", header = FALSE, sep = "\t")
    #message(paste(qval, nrow(bed)))
    bed <- bed[,c(1:3)]
    colnames(bed) <- c("seqnames", "start", "end")
    bed$seqnames <- sub("chr", "", bed$seqnames)
    # write file
    filename <- paste0("data/procdata/tracks/NRF2_", qval, "_sites.bed")
    write.table(
        bed,
        file = filename,
        quote = FALSE,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE
    )
}

# format bed for each qval
qvals <- c("05", "10", "20", "50")
for (qval in qvals) {
    format_bed(qval)
}

# 05 6433
# 10 2555
# 20 679
# 50 467


# -------- from here, run chromvar on cell lines ----------
# pipeline: https://github.com/julianguyn/bca_arche_scoring
# scores for all TF site subsets in NRF2_binding.Zscore.txt
# ---------------------------------------------------------


###########################################################
# Load in data
###########################################################

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# label dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta$sampleid[meta$sampleid %in% dups] <- paste0(meta$sampleid[meta$sampleid %in% dups], " (", meta$tech[meta$sampleid %in% dups], ")")

meta <- meta[meta$type == "cell_line", ]

# chromvar scoring of NRF2 binding
scores <- read.table("data/rawdata/tracks/NRF2_binding.Zscore.txt")

###########################################################
# Format scores (from get_arche_scores())
###########################################################

# standardize sample names of cell lines
colnames(scores) <- sub("^X", "", gsub("\\.(?!$)", "-", colnames(scores), perl = TRUE))
scores <- scores[, colnames(scores) %in% meta$filename]
colnames(scores) <- meta$sampleid[match(colnames(scores), meta$filename)]
scores <- scores[, order(colnames(scores))]

# remove NRF2_50_sites (all NAs)
scores <- scores[-which(rownames(scores) == "NRF2_50_sites"),]

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
    print(Heatmap(scores, cluster_rows = FALSE, name = "NRF2\nBinding\nScore", col = score_pal,
        column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 10), top_annotation = ha))
    dev.off()
}

# plot across all sites
TF_scores_heatmap(scores, "NRF2_all_sites", 3.5)

# plot for individual sites
for (i in 1:nrow(scores)) {
    TF_scores_heatmap(scores[i,], rownames(scores)[i], 2.75)
}