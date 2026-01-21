# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(reshape2)
    library(dplyr)
    library(ComplexHeatmap)
    library(circlize)
    library(genefu)
    library(org.Hs.eg.db)
    library(plyr)
    library(ggplot2)
    library(grid)
    library(gridExtra)
    library(data.table)
})

source("utils/plots/data_exploration.R")
source("utils/score_bca_subtype.R")
source("utils/palettes.R")
source("utils/get_data.R")

###########################################################
# Load in data
###########################################################

# get gene counts from PSets
ubr2 <- get_pset_rna("UBR2")
gray <- get_pset_rna("GRAY")
gcsi <- get_pset_rna("gCSI")
ccle <- get_pset_rna("CCLE")

# get gene counts metadata
meta <- read.table("data/procdata/CCLs/rna/UBR2_RNA_meta.tsv", header = TRUE)

# load in genefu subtyping models
data(pam50.robust)
data(scmgene.robust)
data(scmod2.robust)

###########################################################
# Get overlapping genes
###########################################################

# common genes
genes <- intersect(rownames(ubr2), rownames(ccle))

# subset for overlapping genes
ubr2 <- ubr2[match(genes, rownames(ubr2)),order(colnames(ubr2))]
gray <- gray[match(genes, rownames(gray)),order(colnames(gray))]
gcsi <- gcsi[match(genes, rownames(gcsi)),order(colnames(gcsi))]
ccle <- ccle[match(genes, rownames(ccle)),order(colnames(ccle))]

###########################################################
# Scale RNA-Seq gene counts
###########################################################

# function to scale to 1e6
scale_rna <- function(x) {
  x / sum(x, na.rm = TRUE) * 1e6
}

ubr2 <- apply(ubr2, 2, scale_rna) |> as.data.frame()
gray <- apply(gray, 2, scale_rna) |> as.data.frame()
gcsi <- apply(gcsi, 2, scale_rna) |> as.data.frame()
ccle <- apply(ccle, 2, scale_rna) |> as.data.frame()

###########################################################
# Melt for correlations
###########################################################

# helper function to format for correlations
format_rna <- function(pset) {
    pset$Gene <- rownames(pset)
    pset_m <- reshape2::melt(pset)
    return(pset_m)
}

ubr2_m <- format_rna(ubr2)
gray_m <- format_rna(gray)
gcsi_m <- format_rna(gcsi)
ccle_m <- format_rna(ccle)

###########################################################
# Correlate RNA-seq expression
###########################################################

# helper function to correlate RNA-Seq expression of two psets
corr_pset_rna <- function(pset1, pset2, corr) {

    ccls <- intersect(pset1$variable, pset2$variable)
    pset1 <- pset1[pset1$variable %in% ccls,]
    pset2 <- pset2[pset2$variable %in% ccls,]
    # check that order is the same
    pset1$pairs <- paste0(pset1$Gene, pset1$variable)
    pset2$pairs <- paste0(pset2$Gene, pset2$variable)
    table(pset1$pairs == pset2$pairs)
    corr <- cor(pset1$value, pset2$value,  method = corr, use = "complete.obs")
    return(corr)
}

corr <- "pearson"

p1 <- corr_pset_rna(ubr2_m, gray_m, corr)
p2 <- corr_pset_rna(ubr2_m, gcsi_m, corr)
p3 <- corr_pset_rna(ubr2_m, ccle_m, corr)
p4 <- corr_pset_rna(gray_m, gcsi_m, corr)
p5 <- corr_pset_rna(gray_m, ccle_m, corr)
p6 <- corr_pset_rna(gcsi_m, ccle_m, corr)

###########################################################
# Plot RNA-seq correlation matrix
###########################################################

plot_rna_corr(p1, p2, p3, p4, p5, p6, "Pearson_unlog")

###########################################################
# Compute subtyping model scores
###########################################################

# PAM50
ubr2_pam50 <- score_bca_subtype(ubr2, meta, model = "pam50")
ccle_pam50 <- score_bca_subtype(ccle, meta, model = "pam50")
gcsi_pam50 <- score_bca_subtype(gcsi, meta, model = "pam50")
gray_pam50 <- score_bca_subtype(gray, meta, model = "pam50")

# SCMGENE
ubr2_scmgene <- score_bca_subtype(ubr2, meta, model = "scmgene")

# SCMOD2
ubr2_scmod2 <- score_bca_subtype(ubr2, meta, model = "scmod2")

save(
    ubr2_pam50, ubr2_scmgene, ubr2_scmod2,
    ccle_pam50, gcsi_pam50, gray_pam50,
    file = "data/results/data/3-DataExploration/ccls_subtyping_scores.RData"
)

###########################################################
# Compile and plot across PSets
###########################################################

# get cell line subtype
true_subtype <- get_cell_subtype()

# temp function to compile results across psets
compile_scores <- function(ubr2, gray, gcsi, ccle) {

    # compile subtype results
    toPlot <- rbind.fill(
        as.data.frame(t(ubr2))[rownames(t(ubr2)) == "Subtype",], 
        as.data.frame(t(gray))[rownames(t(gray)) == "Subtype",],
        as.data.frame(t(gcsi))[rownames(t(gcsi)) == "Subtype",],
        as.data.frame(t(ccle))[rownames(t(ccle)) == "Subtype",]
    ) |> t() |> as.data.frame()
    colnames(toPlot) <- c("UBR2", "GRAY", "gCSI", "CCLE")

    toPlot$true_subtype <- true_subtype$subtype[match(rownames(toPlot), true_subtype$sample)]
    return(toPlot)
}

pam50 <- compile_scores(ubr2_pam50, gray_pam50, gcsi_pam50, ccle_pam50)

# plot subtypes across psets
plot_bca_subtype(pam50, "pam50")

###########################################################
# Compile and plot UBR2 across subtyping models
###########################################################

# compile subtype results
toPlot <- rbind.fill(
    as.data.frame(t(ubr2_pam50))[rownames(t(ubr2_pam50)) == "Subtype",],
    as.data.frame(t(ubr2_scmgene))[rownames(t(ubr2_scmgene)) == "Subtype",],
    as.data.frame(t(ubr2_scmod2))[rownames(t(ubr2_scmod2)) == "Subtype",]
) |> t() |> as.data.frame()
colnames(toPlot) <- c("PAM50", "SCMGENE", "SCMOD2")

toPlot$true_subtype <- true_subtype$subtype[match(rownames(toPlot), true_subtype$sample)]

plot_bca_subtype(toPlot, "UBR2")
