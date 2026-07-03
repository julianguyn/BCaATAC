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
    library(matrixStats)
})

source("utils/plots/signatures.R")
source("utils/get_data.R")
source("utils/palettes.R")
source("utils/bca_mutations.R")
source("utils/gene_list.R")
source("utils/plots/ARCHE_scores_heatmap.R")

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

# load in mutation matrix
mut <- read.table("data/results/data/2-MolecularSigAnalysis/TCGA_mutation_matrix.tsv")
meta_mut <- read.csv("metadata/TCGA_mutation_meta.csv")
colnames(mut) <- gsub("\\.", "-", meta_mut$Sample.Name[match(colnames(mut), gsub("-", "\\.", meta_mut$snv_label))])

# load in RNA matrix
rna_df <- as.data.frame(get_tcga_rna())
colnames(rna_df) <- gsub("\\.", "-", colnames(rna_df))

# load in pam50 subtyping
pam50_subtyping <- readRDS("data/procdata/TCGA/pam50_subtyping.rds")

# load in TE zscores
tes <- read.table("data/rawdata/TCGA_TEs/tcga_TE.Zscore.txt")

# load in lehmann classification
# from https://pmc.ncbi.nlm.nih.gov/articles/PMC4911051/#_ad93_
lm <- read_excel("data/rawdata/tcga/pone.0157368.s008.xlsx", sheet = 2, skip = 1, col_names = TRUE)[,1:7] |> as.data.frame()
lm$TCGA_SAMPLE <- sub("-01$", "", lm$TCGA_SAMPLE)
lm <- lm[lm$TCGA_SAMPLE %in% meta$Sample.Name,]

###########################################################
# Handle duplicates and format data
###########################################################

# rename duplicate
mat[mat$variable == "TCGA-A2-A0T4",]$variable <- rep(paste0("TCGA-A2-A0T4-", 1:2), each = 6)
mapping <- unique(data.frame(mat = mat_order, new = mat$variable))
colnames(tfs) <- mapping$new[match(sub("X", "", colnames(tfs)), mapping$mat)]
colnames(tes) <- mapping$new[match(sub("X", "", colnames(tes)), mapping$mat)]

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
pam50 <- format_rna(pam50)

###########################################################
# Format Lehmann data
###########################################################

rownames(lm) <- lm$TCGA_SAMPLE
lm <- as.data.frame(t(lm[,c(6:7)]))

missing <- unique(mat$variable[-which(mat$variable %in% colnames(lm))])
for (sample in missing) {
    lm[[sample]] <- NA
}

lm <- format_rna(lm)

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
    border = border_col
)

###########################################################
# PAM50 heatmap
###########################################################

# make colour palette
cols <- brewer.pal(9, "BuPu")
col_fun <- colorRamp2(
    seq(min(pam50, na.rm = TRUE),
        max(pam50, na.rm = TRUE),
        length.out = 9),
    cols
)

ht2 <- Heatmap(
    pam50,
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
    border = border_col
)

###########################################################
# Lehmann heatmap
###########################################################

ht3 <- Heatmap(
    lm,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "Lehmann\nSubtype",
    column_split = assigned_ARCHE,
    col = lehmann_pal,
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 8),
    row_title = "Lehmann",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10),
    border = border_col
)

###########################################################
# Mutation heatmap
###########################################################

ht4 <- Heatmap(
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
    row_title_gp = gpar(fontsize = 10),
    border = border_col
)

###########################################################
# TEs heatmap
###########################################################

tes <- tes[,match(colnames(toPlot), colnames(tes))] |> as.matrix()
rownames(tes) <- c("", "")

cols <- rev(brewer.pal(9, "PuOr"))
col_fun <- colorRamp2(
    seq(max(tes),
        min(tes),
        length.out = 9),
    cols
)


ht5 <- Heatmap(
    tes,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "Chromvar\nZScore",
    column_split = assigned_ARCHE,
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 8),
    row_title = "TEs",
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10),
    border = border_col
)

###########################################################
# Compiled heatmap
###########################################################

filename <- "data/results/figures/1-Signatures/figure1_heatmap.png"
png(filename, width = 11, height = 8, res = 600, units = "in")
ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5
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

###########################################################
# TFs binding sites heatmap
###########################################################

tfs <- tfs[-which(rownames(tfs) == "IRF8"),]
tfs <- tfs[,match(colnames(toPlot), colnames(tfs))] |> as.matrix()

tfs_norm <- znorm(as.data.frame(tfs))

# subtype annotation
ha1 <- HeatmapAnnotation(
    'ARCHE' = assigned_ARCHE,
    'PAM50' = subtype,
    'gap_spacer' = anno_empty(border = FALSE, height = unit(1, "mm")),
    col = list('ARCHE' = ARCHE_pal, 'PAM50' = subtype_pal),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9)
)


tf_arche <- read.table("data/procdata/TCGA/homer_compiled_TF_genes.txt")
tf_df <- data.frame(
    ARCHE1 = rep(0, nrow(tfs)), ARCHE2 = 0, ARCHE3 = 0,
    ARCHE4 = 0, ARCHE5 = 0, ARCHE6 = 0
)
rownames(tf_df) <- rownames(tfs)
for (gene in rownames(tf_df)) {
    for (arche in paste0("ARCHE", 1:6)) {
        subset_df <- tf_arche[tf_arche$gene == gene,]
        if (arche %in% subset_df$ARCHE) tf_df[gene, arche] <- 1
    }
}


plot_tf_heatmap <- function(toPlot) {

    n <- max(max(toPlot), abs(min(toPlot)))

    # make colour palette
    cols <- rev(brewer.pal(9, "RdBu"))
    col_fun <- colorRamp2(
        seq(5,
            -5,
            length.out = 9),
        cols
    )

    na_val <- "#f8f8f8"
    row_ha <- rowAnnotation(
        ARCHE1 = factor(tf_df$ARCHE1),
        ARCHE2 = tf_df$ARCHE2,
        ARCHE3 = tf_df$ARCHE3,
        ARCHE4 = tf_df$ARCHE4,
        ARCHE5 = tf_df$ARCHE5,
        ARCHE6 = tf_df$ARCHE6,
        col = list(
            'ARCHE1' = c("0" = na_val, "1" = ARCHE_pal[["ARCHE1"]]),
            'ARCHE2' = c("0" = na_val, "1" = ARCHE_pal[["ARCHE2"]]),
            'ARCHE3' = c("0" = na_val, "1" = ARCHE_pal[["ARCHE3"]]),
            'ARCHE4' = c("0" = na_val, "1" = ARCHE_pal[["ARCHE4"]]),
            'ARCHE5' = c("0" = na_val, "1" = ARCHE_pal[["ARCHE5"]]),
            'ARCHE6' = c("0" = na_val, "1" = ARCHE_pal[["ARCHE6"]])
        )
    )


    ht <- Heatmap(
        toPlot,
        cluster_columns = FALSE,
        row_split = 6,
        name = "Chromvar\nZScore",
        column_split = assigned_ARCHE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 8),
        row_names_side = "left",
        column_names_gp = gpar(fontsize = 8),
        row_title = "Top Transcription Factors",
        row_title_side = "left",
        row_title_rot = 90,
        row_title_gp = gpar(fontsize = 10),
        top_annotation = ha1,
        right_annotation = row_ha
    )
    ht

}

plot_tf_heatmap(tfs_norm)
plot_tf_heatmap(tfs)