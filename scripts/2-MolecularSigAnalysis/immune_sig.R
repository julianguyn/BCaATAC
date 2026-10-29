# script to test enrichment of RNA-based immune signatures in ARCHEs

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(GSVA)
    library(GSEABase)
    library(ComplexHeatmap)
    library(circlize)
})

source("utils/get_data.R")
source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# load in metadata
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")

# load in tumour gene counts matrix
t_counts <- get_tcga_rna()

###########################################################
# Get IO signatures
###########################################################

# https://github.com/bhklab/SignatureSets/blob/main/data-raw/signature_information_io.csv
#   TLS_Cabrita: tertiary lymphoid structures
#   IFN_Ayers: IFN-gamma
#   CYT_Rooney, CYT_Davoli: cytolytic activity
#   Teff-IFNG_Fehrenbacher, Teff_McDermott: effector T cells / IFN-gamma
#   CD8_Sade-Feldman: memory-like CD8 T cells
#   T_cells_Danaher_et_al, CD8_T_cells_Danaher_et_al: T cell and CD8 T cell abundance
#   Tcell-inflamed_Ayers: Tumor Inflammation Signature (TIS) genes (Danaher et al. 2018)
#   GeparSixto_Karn: breast cancer immune activation
#   TGFB_Mariathasan: fibroblast TGF-beta / T cell exclusion
#   Exhausted_CD8_Danaher_et_al: exhausted CD8 T cells
#   Macrophage_Danaher, B_cells_Danaher_et_al: macrophage and B cell abundance
#   APM_Wang: antigen processing and presentation

# helper function: mean z-score of signature genes across samples
# genes missing from the matrix or not expressed (zero variance) are dropped
score_signature <- function(expr, genes, name) {
    genes <- intersect(genes, rownames(expr))
    genes <- genes[apply(expr[genes, , drop = FALSE], 1, sd) > 0]
    message(name, ": ", length(genes), " genes used")
    z <- t(scale(t(expr[genes, , drop = FALSE])))
    colMeans(z)
}

sig_files <- list.files("data/rawdata/signatures/IO", pattern = ".csv$", full.names = TRUE)
scores <- data.frame(sample = colnames(t_counts))
for (file in sig_files) {
    sig <- read.csv(file)
    sig <- sig[sig$gene_name != "EIF1AY",] # Y-linked TLS gene, not expressed in this female cohort
    sig <- sig[sig$gene_name != "CD276",] # only negatively weighted TIS gene, excluded from unweighted score
    scores[[sig$signature_name[1]]] <- score_signature(t_counts, sig$gene_name, sig$signature_name[1])
}

###########################################################
# Plot immune signatures by ARCHE
###########################################################

colnames(scores)[colnames(scores) == "Tcell-inflamed_Ayers"] <- "TIS"
arche_meta <- read.csv("data/rawdata/tcga/TCGA_sourcefiles.csv") # ARCHE labels for all tumours (mutation meta only has 70)
scores$ARCHE <- arche_meta$ARCHE[match(scores$sample, arche_meta$Sample.Name)]
scores <- scores[!is.na(scores$ARCHE),]

# mean score of tumours in each ARCHE
toPlot <- aggregate(. ~ ARCHE, data = scores[, -1], FUN = mean)
rownames(toPlot) <- toPlot$ARCHE
toPlot <- t(as.matrix(toPlot[paste0("ARCHE", 1:6), -1]))

# rename signatures
rownames(toPlot) <- sub("_et_al", "", rownames(toPlot))
rownames(toPlot)[rownames(toPlot) == "CYT_Davoli"] <- "CYT-T-NK_Davoli"
rownames(toPlot) <- sub("_[^_]+$", "", rownames(toPlot))

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#B2182B", "#F7F7F7", "#2166AC"))

ha <- HeatmapAnnotation(
    'ARCHE' = colnames(toPlot),
    col = list('ARCHE' = ARCHE_pal),
    simple_anno_size = unit(2, "mm"),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

ht <- Heatmap(
    toPlot,
    name = "Mean\nZ-score",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    rect_gp = gpar(col = "white", lwd = 0.5),
    bottom_annotation = ha,
    heatmap_legend_param = list(
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 7),
        grid_width = unit(3, "mm"),
        grid_height = unit(3, "mm")
    )
)

png("data/results/figures/2-MolecularSigAnalysis/molecularsig/immune_signatures_heatmap.png", width = 2.25, height = 2.75, res = 600, units = "in")
draw(ht)
dev.off()
