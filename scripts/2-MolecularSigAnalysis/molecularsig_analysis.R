# load libraries
suppressPackageStartupMessages({
  library(maftools)
  library(data.table)
  library(GSVA)
  library(GSEABase)
  library(ggplot2)
  library(RColorBrewer)
  library(NMF)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(pheatmap)
  library(dplyr)
  library(qusage)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  library(ggh4x)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
})

source("utils/plots/molecularsig_analysis.R")
source("utils/corr_signatures.R")
source("utils/get_data.R")
source("utils/palettes.R")
source("utils/anno/hallmark_categories.R")
source("utils/anno/sbs_categories.R")

###########################################################
# Load in data
###########################################################

# read in meta data file
meta <- read.csv("metadata/TCGA_mutation_meta.csv")

# load in matrix file from NMF
mat <- get_arche_tcga()
mat <- mat[mat$variable %in% meta$ATAC.Seq.File.Name,]

# get mafs
mafs <- merge_mafs(meta$SNV.File.Name)

# load in tumour gene counts matrix
t_counts <- get_tcga_rna()

# load in ccls gene counts matrix
c_counts <- get_pset_rna("UBR2", gene.symbol = TRUE) |> as.matrix()

# load in gmt files
hallmarks <- getGmt("data/rawdata/gmt/h.all.v2025.1.Hs.symbols.gmt")
myc_targs <- read.gmt("data/rawdata/gmt/All_MYC_Target_Signatures.gmt")

###########################################################
# Extract BCa mutation calls
###########################################################

mafs.tnm <- trinucleotideMatrix(maf = mafs, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
mut_calls <- t(mafs.tnm$nmf_matrix)
colnames(mut_calls) = meta$Sample.Name[match(colnames(mut_calls), meta$snv_label)]
res <- list(signatures = mut_calls)

###########################################################
# Cosine Similarity against Mutational Signatures
###########################################################

# run cosine similarity
og30.cosm <- compareSignatures(nmfRes = res, sig_db = "legacy")$cosine_similarities |> suppressMessages()
v3.cosm <- compareSignatures(nmfRes = res, sig_db = "SBS")$cosine_similarities |> suppressMessages()

# save results
write.csv(og30.cosm, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/og30_cosm.csv", quote = F, row.names = F)
write.csv(v3.cosm, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/v3_cosm.csv", quote = F, row.names = F)

###########################################################
# Single sample gsea on hallmarks gene sets
###########################################################

# ssgsea on tumour counts
t_hm.es <- gsva(ssgseaParam(t_counts, hallmarks), verbose = FALSE)
t_my.es <- gsva(ssgseaParam(t_counts, myc_targs), verbose = FALSE)

# ssgsea on ccls counts
c_hm.es <- gsva(ssgseaParam(c_counts, hallmarks), verbose = FALSE)
c_my.es <- gsva(ssgseaParam(c_counts, myc_targs), verbose = FALSE)

# save results
write.csv(t_hm.es, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/hallmarks_ES.csv", quote = F, row.names = T)
write.csv(t_my.es, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/MYC_ES.csv", quote = F, row.names = T)
write.csv(c_hm.es, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/hallmarks_ccls_ES.csv", quote = F, row.names = T)
write.csv(c_my.es, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/MYC_ccls_ES.csv", quote = F, row.names = T)

###########################################################
# Heatmaps cluster by mutational signatures
###########################################################

plot_pheatmap(og30.cosm, "Original_30")
plot_pheatmap(v3.cosm, "Updated_60")

plot_heatmap_mutsig(og30.cosm, "og30")
plot_heatmap_mutsig(v3.cosm, "v3")

###########################################################
# Correlate each pair of signatures and plot
###########################################################

# format cosine similarity matrices
og30.cosm <- t(og30.cosm) |> as.data.frame()
v3.cosm <- t(v3.cosm) |> as.data.frame()

# correlate signatures
corr_og <- corr_signatures(og30.cosm, "og30")
corr_v3 <- corr_signatures(v3.cosm, "v3")
corr_t_hm <- corr_signatures(t_hm.es, "hm")
corr_t_my <- corr_signatures(t_my.es, "myc")
corr_c_hm <- corr_signatures(c_hm.es, "hm_ccls", "ccls")
corr_c_my <- corr_signatures(c_my.es, "myc_ccls", "ccls")

# plot heatmaps of signature correlations (tumours)
plot_molecularsig_corr(corr_og, "Original 30 COSMIC Signatures", "og30")
plot_molecularsig_corr(corr_v3, "Mutational Signatures", "v3")
plot_molecularsig_corr(corr_t_hm, "Hallmark Gene Sets", "hm")
plot_molecularsig_corr(corr_t_my, "MYC Target Signatures", "myc")
plot_molecularsig_corr(corr_c_hm, "Hallmark Gene Sets", "hm", "ccls")
plot_molecularsig_corr(corr_c_my, "MYC Target Signatures", "myc", "ccls")

# plot box plots
plot_corr_boxplots(corr_v3, corr_t_hm)

###########################################################
# Set up for heatmaps
###########################################################

# make colour palette
cols <- colorRampPalette(c("#39066B", "#B18ED6", "#EBEBEB", "#75BFB2", "#008A8C"))(9)
col_fun <- colorRamp2(seq(-1,1,length.out = 9),cols)

# make row annotation
row_ha <- rowAnnotation(
    ARCHE = paste0("ARCHE", 1:6),
    col = list(
        ARCHE = ARCHE_pal),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

# format toPlot
format_sig_corrs <- function(corr) {
  toPlot <- pivot_wider(
    corr, 
    names_from = Mol.Sig,
    values_from = Corr
  )
  labels <- toPlot$ATAC.Sig
  toPlot$ATAC.Sig <- NULL
  toPlot <- as.data.frame(lapply(toPlot, as.numeric))
  rownames(toPlot) <- labels
  return(toPlot)
}

###########################################################
# Complex heatmap of Hallmarks
###########################################################

# make column annotation
col_ha <- HeatmapAnnotation(
    Annotation = hallmark_categories$category,
    col = list(Annotation = hallmark_pal),
    annotation_name_gp = gpar(fontsize = 9)
)

# format toPlot
toPlot <- format_sig_corrs(corr_t_hm)

# plot
ht <- Heatmap(
    toPlot,
    column_split = 6,
    column_title = "Hallmark Gene Sets",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 10),
    cluster_rows = FALSE,
    name = "Correlation",
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    top_annotation = col_ha,
    right_annotation = row_ha,
    border = border_col,
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (abs(toPlot[i, j]) > 0.5) {
            grid.text("*", x, y, gp = gpar(fontsize = 12, col = "black"))
        }
    }
)

filename <- "data/results/figures/2-MolecularSigAnalysis/molecularsig/figure1_hallmarks_heatmap.png"
png(filename, width = 9, height = 5, res = 600, units = "in")
draw(ht,
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()

###########################################################
# Complex heatmap of Mutational Signatures
###########################################################

# make column annotation
col_ha <- HeatmapAnnotation(
    Annotation = sbs_categories$category,
    col = list(Annotation = sbs_pal),
    annotation_name_gp = gpar(fontsize = 9)
)

# format toPlot
toPlot <- format_sig_corrs(corr_v3)

# plot
ht <- Heatmap(
    toPlot,
    column_split = 6,
    column_title = "Mutational Signatures",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 10),
    cluster_rows = FALSE,
    name = "Correlation",
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    top_annotation = col_ha,
    right_annotation = row_ha,
    border = border_col,
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (abs(toPlot[i, j]) > 0.5) {
            grid.text("*", x, y, gp = gpar(fontsize = 12, col = "black"))
        }
    }
)

filename <- "data/results/figures/2-MolecularSigAnalysis/molecularsig/figure1_cosmic_heatmap.png"
png(filename, width = 9, height = 4.5, res = 600, units = "in")
draw(ht,
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()