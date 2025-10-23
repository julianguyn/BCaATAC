# load libraries
suppressPackageStartupMessages({
    library(readxl)
    library(reshape2)
    library(ggplot2)
    library(ggpubr)
})

source("source/DrugResponse/helper.R")
source("source/DataExploration/helper.R")
source("source/palettes.R")

###########################################################
# Load in TDXd drug response data
###########################################################

# load in TDXd cell line repsonse data
tdxd <- read_excel("DrugResponse/data/TDXd Cell Line Response Data.xlsx", sheet = 1) |>
    as.data.frame()
tdxd <- tdxd[,colnames(tdxd) %in% c("Cell Line", "Average IC50...4", "Average IC50...9")]
colnames(tdxd) <- c("Cell.Line", "Avg.IC50", "Avg.IC50.Treps")

# remove no response values
tdxd$Avg.IC50[tdxd$Avg.IC50 == "#DIV/0!"] <- NA
tdxd$Avg.IC50.Treps[tdxd$Avg.IC50.Treps == "#DIV/0!"] <- NA

###########################################################
# Load in BCa cell line data
###########################################################

# load in BCa cell lines 
samples <- get_cells()

# load in signature scores 
signature_scores <- get_all_scores()

# read in subtype information
true_subtype <- get_cell_subtype()

###########################################################
# Quick cell line mapping
###########################################################

# from map_sen()
for (i in 1:nrow(tdxd)) {
    cell = tdxd$Cell.Line[i]
    if (cell %in% names(mapping_cells)) {tdxd$Cell.Line[i] <- unname(mapping_cells[cell])}
}

###########################################################
# Format TDXd response data 
###########################################################

rownames(tdxd) <- tdxd$Cell.Line
tdxd$Cell.Line <- NULL
tdxd <- t(tdxd) |> as.data.frame()

# keep common samples
tdxd <- tdxd[,which(colnames(tdxd) %in% samples$sample)]
tdxd <- tdxd[,order(colnames(tdxd))]
tdxd_sig <- get_scores(tdxd)

###########################################################
# Stratify HER2 and non-HER2 samples
###########################################################

# get samples
her2 <- true_subtype$sample[true_subtype$subtype == "Her2"]
nonher2 <- true_subtype$sample[true_subtype$subtype != "Her2"]

# subset signatures
her2 <- tdxd_sig[,colnames(tdxd_sig) %in% her2]
nonher2 <- tdxd_sig[,colnames(tdxd_sig) %in% nonher2]

# subset drug response
her2_tdxd <- tdxd[,match(colnames(her2), colnames(tdxd))]
nonher2_tdxd <- tdxd[,match(colnames(nonher2), colnames(tdxd))]

###########################################################
# Load BCa RNA-Seq data
###########################################################

# load in RNA-Seq counts
ubr2 <- load_bca_RNA()

# load in gene metadata
gene_meta <- get_rna_meta()

# keep only ERBB2
keep <- gene_meta$GeneID[gene_meta$Gene.Symbol == "ERBB2"]
erbb2 <- ubr2[rownames(ubr2) == keep,,drop=FALSE] |> as.data.frame()

# keep common cells
common <- intersect(colnames(erbb2), colnames(tdxd))
erbb2 <- erbb2[,match(common, colnames(erbb2))]
erbb2_tdxd <- tdxd[,match(common, colnames(tdxd))]

###########################################################
# Compute Pearson's correlation
###########################################################

tdxd_PC <- computePC(tdxd_sig, tdxd, "TDXd")
#her2_PC <- computePC(her2, her2_tdxd, "TDXd")              # only 3 HER2 CCLs
#nonher2_PC <- computePC(nonher2, nonher2_tdxd, "TDXd")
erbb2_PC <- computePC(erbb2, erbb2_tdxd, "TDXd")

###########################################################
# Combine data for plotting
###########################################################

all_c <- combine_data(tdxd_sig, tdxd)
her_c <- combine_data(her2, her2_tdxd)
non_c <- combine_data(nonher2, nonher2_tdxd)
erbb2_c <- combine_data(erbb2, erbb2_tdxd)

###########################################################
# Plot correlation of ARCHE and TDXd response
###########################################################

plot_tdxd_all(all_c, tdxd_PC, "Avg.IC50", "All Cells")
plot_tdxd_all(all_c, tdxd_PC, "Avg.IC50.Treps", "All Cells")
plot_tdxd_all(her_c, tdxd_PC, "Avg.IC50", "HER2 Cells")
plot_tdxd_all(her_c, tdxd_PC, "Avg.IC50.Treps", "HER2 Cells")
plot_tdxd_all(non_c, tdxd_PC, "Avg.IC50", "Non-HER2 Cells")
plot_tdxd_all(non_c, tdxd_PC, "Avg.IC50.Treps", "Non-Her2 Cells")

plot_tdxd_corr(erbb2_c, erbb2_PC, "ENSG00000141736.14", "Avg.IC50", save = TRUE)
plot_tdxd_corr(erbb2_c, erbb2_PC, "ENSG00000141736.14", "Avg.IC50.Treps", save = TRUE)