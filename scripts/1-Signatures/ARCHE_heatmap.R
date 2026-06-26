# main figure to plot heatmap of ARCHE and other molecular features

suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
    library(tidyverse)
    library(patchwork)
    library(RColorBrewer)
    library(reshape2)
    library(ComplexHeatmap)
    library(circlize)
})

source("utils/plots/signatures.R")
source("utils/get_data.R")
source("utils/palettes.R")
source("utils/bca_mutations.R")
source("utils/gene_list.R")

set.seed(123)

###########################################################
# Load in data
###########################################################

# load in metadata
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
meta$Sample.Name <- gsub("\\.", "-", meta$Sample.Name)

# get signature scores from NMF
mat <- get_arche_tcga()
mat$variable <- meta$Sample.Name[match(mat$variable, meta$ATAC.Seq.File.Name)]

# load in mutation matrix
mut <- read.table("data/results/data/2-MolecularSigAnalysis/TCGA_mutation_matrix.tsv")
meta <- read.csv("metadata/TCGA_mutation_meta.csv")
colnames(mut) <- gsub("\\.", "-", meta$Sample.Name[match(colnames(mut), gsub("-", "\\.", meta$snv_label))])

# load in RNA matrix
rna_df <- as.data.frame(get_tcga_rna())
colnames(rna_df) <- gsub("\\.", "-", colnames(rna_df))

# load in pam50 subtyping
pam50_subtyping <- readRDS("data/procdata/TCGA/pam50_subtyping.rds")

# load in correlated genes
#gene_corr <- readRDS("data/results/data/2-MolecularSigAnalysis/gene_correlations_arches.rds")

###########################################################
# Get DEGS
###########################################################

genes <- c()
for (arche in paste0("ARCHE", c(1,4,6,2,5,3))) {
    deg <- read.csv(paste0("data/results/data/2-MolecularSigAnalysis/DEG/ARCHE_", arche, "_vs_Other.csv"))
    deg <- deg[which(deg$padj < 0.05),]
    deg <- deg[order(abs(deg$log2FoldChange), decreasing = TRUE),]
    #print(nrow(deg[abs(deg$log2FoldChange) > 7,]))
    genes <- c(genes, deg$X[abs(deg$log2FoldChange) > 6])
    genes <- unique(genes)
}

###########################################################
# Handle duplicates and format data
###########################################################

# rename duplicate
mat[mat$variable == "TCGA-A2-A0T4",]$variable <- rep(paste0("TCGA-A2-A0T4-", 1:2), each = 6)

# format toPlot
toPlot <- mat %>% select(ARCHE, value, variable)
toPlot <- toPlot %>%
    pivot_wider(
        names_from = c(variable),
        values_from = value
    ) %>%
    column_to_rownames(var = "ARCHE")

###########################################################
# Get variables
###########################################################

assigned_ARCHE <- unique(mat[,c(2,4)])$signature_assign
subtype <- unique(mat[,c(2,5)])$subtype

###########################################################
# Format mutation data
###########################################################

missing <- unique(mat$variable[-which(mat$variable %in% colnames(mut))])
for (sample in missing) {
    mut[[sample]] <- NA
}
mut$'TCGA-A2-A0T4-1' <- mut$'TCGA-A2-A0T4-2' <- mut[,colnames(mut) == "TCGA-A2-A0T4"]
mut$'TCGA-A2-A0T4' <- NULL
mut <- mut[rownames(mut) %in% bca_mutations,]

mut <- mut[,match(colnames(toPlot), colnames(mut))] |> as.matrix()
mut_bin <- ifelse(is.na(mut), NA, ifelse(mut >= 1, 1, 0))

###########################################################
# Format RNA data
###########################################################

format_rna <- function(rna) {
    rna$'TCGA-A2-A0T4-1' <- rna$'TCGA-A2-A0T4-2' <- rna[,colnames(rna) == "TCGA-A2-A0T4"]
    rna$'TCGA-A2-A0T4' <- NULL

    rna <- rna[,match(colnames(toPlot), colnames(rna))] |> as.matrix()
    return(rna)
}

pam50 <- as.data.frame(t(pam50_subtyping$subtype.proba))
rna <- format_rna(pam50)

###########################################################
# ARCHE heatmap
###########################################################

# make colour palette
cols <- brewer.pal(9, "Blues")
col_fun <- colorRamp2(
    seq(min(toPlot, na.rm = TRUE),
        max(toPlot, na.rm = TRUE),
        length.out = 9),
    cols
)

# subtype annotation
ha1 <- HeatmapAnnotation(
    'ARCHE' = assigned_ARCHE,
    'PAM50' = subtype,
    'gap_spacer' = anno_empty(border = FALSE, height = unit(1, "mm")),
    col = list('ARCHE' = ARCHE_pal, 'PAM50' = subtype_pal),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9)
)

# heatmap
ht1 <- Heatmap(
    toPlot,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "NMF\nWeight",
    column_split = assigned_ARCHE,
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 8),
    top_annotation = ha1,
    row_title = "ARCHEs",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10),
    border = TRUE
)

###########################################################
# PAM50 heatmap
###########################################################

# make colour palette
cols <- brewer.pal(9, "BuPu")
col_fun <- colorRamp2(
    seq(min(rna, na.rm = TRUE),
        max(rna, na.rm = TRUE),
        length.out = 9),
    cols
)

ht2 <- Heatmap(
    rna,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "Subtyping\nScore",
    column_split = assigned_ARCHE,
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 8),
    row_title = "PAM50",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10),
    border = TRUE
)

###########################################################
# Mutation heatmap
###########################################################

ht3 <- Heatmap(
    mut_bin,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "Mutation\nStatus",
    column_split = assigned_ARCHE,
    col = c("0" = "#F1F1F1", "1" = random_blue),
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 8),
    row_title = "Mutations",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10)
)

###########################################################
# Compiled heatmap
###########################################################

filename <- "data/results/figures/1-Signatures/figure1_heatmap.png"
png(filename, width = 11, height = 8, res = 600, units = "in")
ht1 %v% ht2 %v% ht3
dev.off()

###########################################################
# PAM50 gene heatmap
###########################################################

rna_df <- rna_df[pam50_genes,]
rna_df <- format_rna(rna_df)

# make colour palette
cols <- colorRampPalette(c("#7F2929", "#BA5050", "#D4D7F8", "#80AAEE", "#3878DF"))(9)
col_fun <- colorRamp2(
    seq(0,
        20,
        length.out = 9),
    cols
)

# subtype annotation
ha1 <- HeatmapAnnotation(
    'ARCHE' = assigned_ARCHE,
    'PAM50' = subtype,
    'gap_spacer' = anno_empty(border = FALSE, height = unit(1, "mm")),
    col = list('ARCHE' = ARCHE_pal, 'PAM50' = subtype_pal),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9)
)

ht4 <- Heatmap(
    rna_df,
    row_split = 8,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    name = "Gene\nExpression",
    column_split = assigned_ARCHE,
    col = col_fun,
    top_annotation = ha1,
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 8),
    row_title = "PAM50 Genes",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10)
)

filename <- "data/results/figures/1-Signatures/figure1_pam50_gene_heatmap.png"
png(filename, width = 10, height = 8, res = 600, units = "in")
ht4
dev.off()

###########################################################
# Claudin heatmap
###########################################################

rna_cl_u <- format_rna(as.data.frame(rna_df[claudin_low_up,]))
rna_cl_d <- format_rna(as.data.frame(rna_df[claudin_low_down,]))

# make colour palette
cols <- brewer.pal(9, "Purples")
col_fun <- colorRamp2(
    seq(min(toPlot, na.rm = TRUE),
        max(toPlot, na.rm = TRUE),
        length.out = 9),
    cols
)

# subtype annotation
ha1 <- HeatmapAnnotation(
    'ARCHE' = assigned_ARCHE,
    'PAM50' = subtype,
    col = list('ARCHE' = ARCHE_pal, 'PAM50' = subtype_pal),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9)
)

ht1 <- Heatmap(
    rna_cl_u,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    name = "Gene\nExpression",
    column_split = assigned_ARCHE,
    col = col_fun,
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_names_gp = gpar(fontsize = 8),
    row_title = "CL (upregulated)",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    top_annotation = ha1
)

ht2 <- Heatmap(
    rna_cl_d,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    name = "Gene\nExpression",
    column_split = assigned_ARCHE,
    col = col_fun,
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_names_gp = gpar(fontsize = 8),
    row_title = "CL (downregulated)",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 12, fontface = "bold")
)

filename <- "data/results/figures/1-Signatures/suppfigure1_cl_heatmap.png"
png(filename, width = 9, height = 7, res = 600, units = "in")
draw(ht1 %v% ht2, merge_legends = TRUE)
dev.off()