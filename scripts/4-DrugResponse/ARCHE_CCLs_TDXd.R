###########################################################
# Format TDXd response data
###########################################################

# get tdxd response data
tdxd <- get_tdxd()

# map cell lines names
rownames(tdxd) <- map_cells(tdxd$Cell.Line)
tdxd$Cell.Line <- NULL
tdxd <- t(tdxd) |> as.data.frame()

# keep common samples
tdxd <- tdxd[,which(colnames(tdxd) %in% c_meta$sampleid)]
tdxd_sen <- tdxd[,order(colnames(tdxd))]

###########################################################
# Load in gene count data (for tdxd analysis)
###########################################################

# get BCa gene counts
ubr2 <- get_pset_rna("UBR2")

# get gene counts metadata
gene_meta <- read.table("data/procdata/CCLs/rna/UBR2_RNA_meta.tsv", header = TRUE)

# keep only ERBB2
keep <- gene_meta$Ensembl[gene_meta$Gene.Symbol == "ERBB2"]
erbb2 <- ubr2[rownames(ubr2) == keep,,drop=FALSE] |> as.data.frame()

# keep common cells
common <- intersect(colnames(erbb2), colnames(tdxd_sen))
erbb2 <- erbb2[,match(common, colnames(erbb2))]
erbb2_tdxd <- tdxd_sen[,match(common, colnames(tdxd_sen))]

###########################################################
# Stratify HER2 and non-HER2 samples (for tdxd analysis)
###########################################################

# get samples
her2 <- c_meta$sampleid[c_meta$subtype == "Her2"]
nonher2 <- c_meta$sampleid[c_meta$subtype != "Her2"]

# subset signatures
her2 <- tdxd_sig[,colnames(tdxd_sig) %in% her2]
nonher2 <- tdxd_sig[,colnames(tdxd_sig) %in% nonher2]

# subset drug response
her2_tdxd <- tdxd_sen[,match(colnames(her2), colnames(tdxd_sen))]
nonher2_tdxd <- tdxd_sen[,match(colnames(nonher2), colnames(tdxd_sen))]

# todo:: implement tdxd associations with new arche scores
tdxd_PC <- compute_pc(tdxd_sig, tdxd_sen, "TDXd")
erbb2_PC <- compute_pc(erbb2, erbb2_tdxd, "TDXd")

###########################################################
# Combine TDXd data for plotting
###########################################################

# helper function to combine drug response and arche scores
combine_data <- function(signature_scores, tdxd) {
    combined <- t(rbind(tdxd, signature_scores)) |> as.data.frame()
    samples <- rownames(combined)
    combined <- sapply(combined, as.numeric) |> as.data.frame()
    rownames(combined) <- samples
    combined$subtype <- true_subtype$subtype[match(samples, true_subtype$sample)]
    return(combined)
}

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

plot_tdxd_corr(erbb2_c, erbb2_PC, "ENSG00000141736", "Avg.IC50", save = TRUE)
plot_tdxd_corr(erbb2_c, erbb2_PC, "ENSG00000141736", "Avg.IC50.Treps", save = TRUE)