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

source("source/DataExploration/helper.R")
source("source/DataExploration/plots.R")
source("source/DrugResponse/helper.R")
source("source/palettes.R")
source("source/get_data.R")

###########################################################
# Load in data
###########################################################

# get gene counts from PSets
ubr2 <- get_pset_rna("UBR2")
gray <- get_pset_rna("GRAY")
gcsi <- get_pset_rna("gCSI")
ccle <- get_pset_rna("CCLE")

# get gene counts metadata
meta <- read.table("Preprocessing/procdata/CCLs/UBR2_RNA_meta.tsv", header = TRUE)

# load in genefu subtyping models
data(pam50.robust)
data(scmgene.robust)
data(scmod2.robust)

###########################################################
# Melt for correlations
###########################################################

# common genes
genes <- intersect(rownames(ubr2), rownames(ccle))

# melt
ubr2_m <- reshape2::melt(ubr2[match(genes, rownames(ubr2)),order(colnames(ubr2))])
gray_m <- reshape2::melt(gray[match(genes, rownames(gray)),order(colnames(gray))])
gcsi_m <- reshape2::melt(gcsi[match(genes, rownames(gcsi)),order(colnames(gcsi))])
ccle_m <- reshape2::melt(ccle[match(genes, rownames(ccle)),order(colnames(ccle))])


###########################################################
# Correlate RNA-seq expression
###########################################################

p1 <- corr_pset_rna(ubr2_m, gray_m)
p2 <- corr_pset_rna(ubr2_m, gcsi_m)
p3 <- corr_pset_rna(ubr2_m, ccle_m)
p4 <- corr_pset_rna(gray_m, gcsi_m)
p5 <- corr_pset_rna(gray_m, ccle_m)
p6 <- corr_pset_rna(gcsi_m, ccle_m)

###########################################################
# Plot RNA-seq correlation matrix
###########################################################

plot_rna_corr(p1, p2, p3, p4, p5, p6)

###########################################################
# Compute subtyping model scores
###########################################################

# PAM50
ubr2_pam50 <- score_bcasubtype(ubr2, meta, model = "pam50")
ccle_pam50 <- score_bcasubtype(ccle, meta, model = "pam50")
gcsi_pam50 <- score_bcasubtype(gcsi, meta, model = "pam50")
gray_pam50 <- score_bcasubtype(gray, meta, model = "pam50")

# SCMGENE
ubr2_scmgene <- score_bcasubtype(ubr2, meta, model = "scmgene")

# SCMOD2
ubr2_scmod2 <- score_bcasubtype(ubr2, meta, model = "scmod2")

save(
    ubr2_pam50, ubr2_scmgene, ubr2_scmod2,
    ccle_pam50, gcsi_pam50, gray_pam50, 
    file = "DataExploration/results/data/subtyping_scores.RData"
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
plot_bcasubtype(pam50, "pam50")

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

plot_bcasubtype(toPlot, "UBR2")
