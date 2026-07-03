# check for stemness TE signatures

# analyze individual groups of TEs

suppressPackageStartupMessages({
    library(readxl)
    library(tidyverse)
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    library(reshape2)
})

source("utils/get_data.R")
source("utils/palettes.R")

set.seed(123)

## --------------------------------------------------------
# Preprocessing for H4H

# load in stemness signatures
te_tp <- read_excel("data/rawdata/TEs/cd-23-0331_supplementary_table_s2_suppst2.xlsx", sheet = 1, col_names = TRUE, skip = 1)
te_pt <- read_excel("data/rawdata/TEs/cd-23-0331_supplementary_table_s3_suppst3.xlsx", sheet = 1, col_names = TRUE, skip = 1)

te_tp <- te_tp$'TE family'[te_tp$Enrichment == "Mammary"] #56
te_pt <- te_pt$'TE family'[te_pt$Enrichment == "PSCs_Mammary"] #194

save(te_tp, te_pt, file = "data/procdata/TCGA/TE_stemness_signature.RData")

# from here, run 0-Preprocessing/process_TEs.R

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
tes <- read.table("data/rawdata/TCGA_families/tcga_TE_families_all.Zscore.txt")

# load in stemness signatures
te_tp <- read_excel("data/rawdata/TEs/cd-23-0331_supplementary_table_s2_suppst2.xlsx", sheet = 1, col_names = TRUE, skip = 1)
te_pt <- read_excel("data/rawdata/TEs/cd-23-0331_supplementary_table_s3_suppst3.xlsx", sheet = 1, col_names = TRUE, skip = 1)

###########################################################
# Get TEs
###########################################################

te_tp <- te_tp$'TE family'[te_tp$Enrichment == "Mammary"] #56
te_pt <- te_pt$'TE family'[te_pt$Enrichment == "PSCs_Mammary"] #194

save(te_tp, te_pt, file = "data/procdata/TE_stemness_signature.RData")

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

create_TE_toPlot <- function(tes) {

    toPlot <- data.frame(matrix(nrow = nrow(tes), ncol = 6))
    rownames(toPlot) <- rownames(tes)
    colnames(toPlot) <- paste0("ARCHE", 1:6)

    for (arche in colnames(toPlot)) {
        samples <- tcga_anno$variable[tcga_anno$signature_assign == arche]
        tt <- tes[,which(colnames(tes) %in% samples)]
        toPlot[[arche]] <- rowMeans(tt)
    }
    return(toPlot)
}

te_tp <- create_TE_toPlot(tes[te_tp,])
te_pt <- create_TE_toPlot(tes[te_pt,])

###########################################################
# Colour palettes
###########################################################

cols <- colorRampPalette(c("#4B8F8C", "#A2E4D0","#AC71A7", "#5E4352"))(9)
n <- 2
col_fun <- colorRamp2(seq(n, -n, length.out = 9), cols)

###########################################################
# Plot heatmaps
###########################################################

ht1 <- Heatmap(
    te_tp,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    name = "Capped\nZscore",
    col = col_fun,
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_names_gp = gpar(fontsize = 8),
    row_title = "Tissue-state\nEnriched TEs",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 0,
    column_names_centered = TRUE
)

ht2 <- Heatmap(
    te_pt,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    name = "Capped\nZscore",
    col = col_fun,
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_names_gp = gpar(fontsize = 8),
    row_title = "Pluripotent Stem Cell\nEnriched TEs",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 0,
    column_names_centered = TRUE,
)

filename <- "data/results/figures/1-Signatures/suppfigure1_te_stemness_heatmap.png"
png(filename, width = 5, height = 7, res = 600, units = "in")
draw(ht1 %v% ht2, merge_legends = TRUE)
dev.off()

###########################################################
# Plot violin plots
###########################################################

plot_violin <- function(te) {
    te$TE <- rownames(te)
    toPlot <- reshape2::melt(te)

    ggplot(toPlot, aes(x = variable, y = value, fill = variable)) + 
        geom_violin() +
        geom_boxplot(width = 0.3) +
        geom_jitter(
            position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2),
            alpha = 0.2, size = 1
        ) +
        scale_fill_manual(values = ARCHE_pal) +
        theme_bw()
}