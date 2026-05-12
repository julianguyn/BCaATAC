# umap - run on H4H

library(data.table)

INDIR <- "/cluster/projects/bhklab/projects/BCaATAC/BCa_ARCHE_Scoring/data/results/"

###########################################################
# Load in data
###########################################################

# helper function to read in binary matrix
load_mat <- function(path) {
    mat <- readRDS(paste0(INDIR, path))
    setDT(mat)
    return(mat)
}

# load in cell lines
cell <- load_mat("cell_lines/50k_TCGAbg/matrix/cells_50k.consensus.Binarymat.rds")
tcga_cells <- load_mat("TCGA/cell_50k/matrix/TCGA_50k_cell.consensus.Binarymat.rds")

# load in pdxs
pdx <- load_mat("PDXs/50k_TCGAbg/matrix/PDXs_50k.consensus.Binarymat.rds")
tcga_pdxs <- load_mat("TCGA/pdx_50k/matrix/TCGA_50k_pdx.consensus.Binarymat.rds")

###########################################################
# Merge matrices and format for pca
###########################################################

# helper function to format dataframes for pca
format_df <- function(sample, tcga) {

  # merged on peaks
  merged <- merge(
    cell,
    tcga,
    by = c("seqnames", "start", "end"),
    all = TRUE
  )

  # replace NA values with 0
  for (j in names(merged)) {
    set(merged,
        which(is.na(merged[[j]])),
        j,
        0)
  }

  # format for PCA
  merged <- as.data.frame(merged)
  rownames(merged) <- paste0(merged$seqnames, ":", merged$start, ":", merged$end)
  merged <- merged[,-c(1:3)]
  merged <- t(merged)

  return(merged)
}

merged_cells <- format_df(cells, tcga_cells)
merged_pdxs <- format_df(pdxs, tcga_pdxs)

###########################################################
# Run PCA
###########################################################

# helper function to run pca
run_pca <- function(merged, outfile) {

  # run PCA
  pca_res <- prcomp(merged)
  samples <- rownames(merged)
  print(head(samples))

  # get variances
  variance_pca <- pca_res$sdev^2
  variance_pca / sum(variance_pca) -> prop_var
  print(head(prop_var))

  pca_res <- pca_res$x[,c(1:2)]
  print(head(pca_res))

  save(pca_res, samples, prop_var, file = outfile)

}

run_pca(merged_cells, "cell_tcga_pca.RData")
run_pca(merged_pdxs, "pdx_tcga_pca.RData")

###########################################################
# Plot PCA
###########################################################

load("data/procdata/scoring/cell_tcga_pca.RData")

plot(prop_var,
     xlab = "Principal Components",
     ylab = "Proportion of Variance Explained",
     ylim = c(0, 1))
