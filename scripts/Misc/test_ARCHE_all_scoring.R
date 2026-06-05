# test ARCHE scores with all samples + TCGA (cell + pdx + tcga)

# load libraries
suppressPackageStartupMessages({
    library(matrixStats)
    library(data.table)
    library(ComplexHeatmap)
    library(circlize)
})

set.seed(101)

source("utils/palettes.R")
source("utils/plots/ARCHE_scores_heatmap.R")

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
# Try removing low deviation samples
###########################################################
### - removed

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
        Subtype = meta$subtype[match(colnames(df), meta$sampleid)],
        Type = meta$type[match(colnames(df), meta$sampleid)],
        SumDevs = colSums(abs(df)),
        col = list(Subtype = subtype_pal, Type = sample_type_pal, SumDevs = count_pal))

    filename <- paste0("data/results/figures/Misc/all_scoring/", label, "_ARCHE_scores.png")
    cat("-----Saving plot to", filename, "\n")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(
        Heatmap(toPlot, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 3),
            row_names_gp = gpar(fontsize = 10), top_annotation = ha)
    )
    dev.off()
}

# plot heatmap of zscores
plot_ARCHE_scores_heatmap_counts(zsc, "zcores", meta)
plot_ARCHE_scores_heatmap_counts(zsc, "znorm_zscores", meta, znorm = TRUE)

# heatmap of deviations
plot_ARCHE_scores_heatmap_counts(dev, "deviations", meta)
plot_ARCHE_scores_heatmap_counts(dev, "znorm_deviations", meta, znorm = TRUE)
