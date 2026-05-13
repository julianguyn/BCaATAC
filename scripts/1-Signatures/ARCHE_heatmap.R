# main figure to plot heatmap of ARCHE and other molecular features

suppressPackageStartupMessages({
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
    'PAM50\nSubtype' = subtype,
    col = list('PAM50\nSubtype' = subtype_pal),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9)
)

# heatmap
ht1 <- Heatmap(
    toPlot,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "ARCHE\nScore",
    column_split = assigned_ARCHE,
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    top_annotation = ha1
)

###########################################################
# Mutation heatmap
###########################################################

ht2 <- Heatmap(
    mut_bin,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "Mutation\nStatus",
    column_split = assigned_ARCHE,
    col = c("0" = "#F1F1F1", "1" = random_blue),
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_names_gp = gpar(fontsize = 8)
)

###########################################################
# Compiled heatmap
###########################################################

filename <- "data/results/figures/1-Signatures/"
png(filename, width = 11, height = 6, res = 600, units = "in")
ht1 %v% ht2
dev.off()
