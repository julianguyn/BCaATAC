# script to format tracks from Chip-Atlas for chromVar

suppressPackageStartupMessages({
    library(data.table)
    library(matrixStats)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
})

###########################################################
# Format TF binding sites bed file
###########################################################

# load in and format file
bed <- fread("data/rawdata/tracks/Oth.Brs.05.NFE2L2.AllCell.bed", skip = "track", header = FALSE, sep = "\t")
bed <- bed[,c(1:3)]
colnames(bed) <- c("seqnames", "start", "end")
bed$seqnames <- sub("chr", "", bed$seqnames)

# write file
filename <- paste0("data/rawdata/tracks/NRF2_TF_sites.bed")
write.table(bed, file = filename, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

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

###########################################################
# Plot output of chromvar (from plots/ARCHE_scores_heatmap.R)
###########################################################

rownames(scores) <- "NRF2 TF Score"

# set colours for plotting
lim <- max(c(abs(min(scores)), max(scores)))
score_pal <- colorRamp2(seq(-lim, lim, length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))

ha <- HeatmapAnnotation(
    Subtype = meta[match(colnames(scores), meta$sampleid),]$subtype,
    col = list(Subtype = subtype_pal)
)

filename <- paste0("data/results/figures/Misc/NRF2_TF_scores.png")
png(filename, width = 10, height = 2.75, res = 600, units = "in")
Heatmap(scores, cluster_rows = FALSE, name = "NRF2\nBinding\nScore", col = score_pal,
    column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
    row_names_gp = gpar(fontsize = 10), top_annotation = ha)
dev.off()
