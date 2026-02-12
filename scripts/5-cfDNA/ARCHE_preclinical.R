# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
})

source("utils/score_arche_cfDNA.R")
source("utils/plots/cfdna.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/palettes.R")
source("utils/mappings.R")

###########################################################
# Get ARCHE scores
###########################################################

dir <- "preclinical-ARCHE"
scores <- score_arche_cfDNA(dir)
scores$ARCHE <- sub("sig", "ARCHE", scores$ARCHE)

###########################################################
# Map sample names
###########################################################

scores$Sample <- factor(
  scores$Sample, 
  levels = c("CAMA1_30", "MCF7_merged", "BPTO95", "CAMA1_mouse_ctDNA_merged"),
  labels = c("CCL (CAMA-1)", "CCL (MCF-7)", "Organoid", "Xenograft")
)

###########################################################
# Plot ARCHE scores
###########################################################

plot_ARCHE_score_preclinical(scores)

###########################################################
# Prepare ARCHE scores of matched ATAC-Seq
###########################################################

# look in ARCHE scores manually
atac_scores <- read.table("data/rawdata/cfDNA-ATAC/chromvar/cfDNA_ATAC.Zscore.txt")
rownames(atac_scores) <- paste0("ARCHE", 1:6)
colnames(atac_scores) <- sub("ATAC_", "", sub("_ATAC", "", sub("_S.*", "", colnames(atac_scores))))

# create metadata
reps <- c("Rep1", "Rep2")
meta <- data.frame(
  SampleID = colnames(atac_scores),
  Sample = sub("_.*", "", sub("RES.*", "", sub("PAR.*", "", colnames(atac_scores)))),
  Model = c(rep("CCL", 12), "PDO", "PDX", "PDX", "PDO", rep("CCL", 8)),
  Res = "Parental",
  Rep = c(rep(reps, 6), "NA", reps, "NA", rep(reps, 4))
)
meta$Res[grep("RES", colnames(atac_scores))] <- "Resistant"

###########################################################
# Manually plot heatmap
# code from: plot_ARCHE_scores_heatmap()
###########################################################

# normalize scores
df <- znorm(atac_scores)

# set colours for plotting
lim <- max(c(abs(min(df)), max(df)))
score_pal = colorRamp2(seq(-lim, lim, length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))

ha <- HeatmapAnnotation(
  Sample = meta[match(colnames(df), meta$SampleID),]$Sample,
  Model = meta[match(colnames(df), meta$SampleID),]$Model,
  Res = meta[match(colnames(df), meta$SampleID),]$Res,
  Rep = meta[match(colnames(df), meta$SampleID),]$Rep
)

filename <- paste0("data/results/figures/5-cfDNA/ARCHEheatmaps/sampleOnly_ARCHE_scores_norm.png")
png(filename, width = 8, height = 4, res = 600, units = "in")
print(
    Heatmap(df, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
        column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 10), top_annotation = ha)
)
dev.off()