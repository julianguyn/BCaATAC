# analyze individual groups of TEs

suppressPackageStartupMessages({
    library(tidyverse)
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
})

source("utils/get_data.R")
source("utils/palettes.R")

set.seed(123)

###########################################################
# Load in data
###########################################################

# load in metadata
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
meta$Sample.Name <- gsub("\\.", "-", meta$Sample.Name)

# get signature scores from NMF
mat <- get_arche_tcga()
mat_order <- mat$variable
mat$variable <- meta$Sample.Name[match(mat$variable, meta$ATAC.Seq.File.Name)]

# load in HERV_LTRs zscores
tes <- read.table("data/rawdata/TCGA_HERV_LTRs/tcga_HERV_LTRs.Zscore.txt")

###########################################################
# Handle duplicates and format data
###########################################################

# rename duplicate
mat[mat$variable == "TCGA-A2-A0T4",]$variable <- rep(paste0("TCGA-A2-A0T4-", 1:2), each = 6)
mapping <- unique(data.frame(mat = mat_order, new = mat$variable))
colnames(tes) <- mapping$new[match(sub("X", "", colnames(tes)), mapping$mat)]

###########################################################
# Get variables
###########################################################

tcga_anno <- unique(mat[,c(2,4,5)])

###########################################################
# Create toPlot
###########################################################

toPlot <- data.frame(matrix(nrow = nrow(tes), ncol = 6))
rownames(toPlot) <- rownames(tes)
colnames(toPlot) <- paste0("ARCHE", 1:6)

for (arche in colnames(toPlot)) {
    samples <- tcga_anno$variable[tcga_anno$signature_assign == arche]
    tt <- tes[,which(colnames(tes) %in% samples)]
    toPlot[[arche]] <- rowMeans(tt)
}

###########################################################
# Plot
###########################################################

# make colour palette
cols <- rev(brewer.pal(9, "RdBu"))
col_fun <- colorRamp2(
    seq(8,
        -8,
        length.out = 9),
    cols
)

# ARCHE annotation
ha1 <- HeatmapAnnotation(
    'ARCHE' = paste0("ARCHE", 1:6),
    'gap_spacer' = anno_empty(border = FALSE, height = unit(1, "mm")),
    col = list('ARCHE' = ARCHE_pal),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9)
)

ht <- Heatmap(
    toPlot,
    cluster_columns = FALSE,
    name = "Chromvar\nZScore",
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 0,
    column_names_centered = TRUE,
    row_title = "TE Family",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10),
    top_annotation = ha1
)
ht
h