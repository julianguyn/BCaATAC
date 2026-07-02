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

# load in zscores
herv_ltrs <- read.table("data/rawdata/TCGA_families/tcga_HERV_LTRs.Zscore.txt")
tes <- read.table("data/rawdata/TCGA_families/tcga_TE_families.Zscore.txt")
counts <- readRDS("data/rawdata/TCGA_families/count_regions.rds")

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

toPlot <- toPlot[rownames(toPlot) %in% counts$Family[counts$Regions > 5000],]

###########################################################
# Colour palettes
###########################################################

cols <- colorRampPalette(c("#4B8F8C", "#A2E4D0","#AC71A7", "#5E4352"))(9)
n <- 5
col_fun <- colorRamp2(seq(n, -n, length.out = 9), cols)

herv_pal <- c("Yes" = "#0F5257", "No" = "grey85")

log2count_col <- colorRamp2(
  seq(min(5, na.rm = TRUE), max(log2(counts$Regions), na.rm = TRUE), length.out = 9),
  brewer.pal(9, "BuPu")
)

###########################################################
# Create row anno
###########################################################

counts <- counts[match(rownames(toPlot), counts$Family),]
counts$HERV_LTRs <- ifelse(counts$Family %in% rownames(herv_ltrs), "Yes", "No")

row_ha <- rowAnnotation(
    Log2Count = log2(counts$Regions), HERV_LTRs = counts$HERV_LTRs,
    col = list(
        Log2Count = log2count_col,
        HERV_LTRs = herv_pal
    ),
    annotation_name_gp = gpar(fontsize = 8)
)


###########################################################
# Plot
###########################################################

col_ha <- HeatmapAnnotation(
    'ARCHE' = paste0("ARCHE", 1:6),
    col = list('ARCHE' = ARCHE_pal),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

ht <- Heatmap(
    toPlot,
    cluster_columns = FALSE,
    name = "Capped\nChromvar\nZScore",
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
    bottom_annotation = col_ha,
    right_annotation = row_ha
)
filename <- "data/results/figures/1-Signatures/TE_heatmap.png"
png(filename, width = 6.5, height = 5, res = 600, units = "in")
ht
dev.off()
