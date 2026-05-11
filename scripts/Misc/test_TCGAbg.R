# umap - run on H4H

library(data.table)

INDIR <- "/cluster/projects/bhklab/projects/BCaATAC/BCa_ARCHE_Scoring/data/results/"

# helper function to read in binary matrix
load_mat <- function(path) {
    mat <- readRDS(paste0(INDIR, path))
    mat <- mat[1:50,1:50]
    setDT(mat)
    return(mat)
}

# load in objects
cell <- load_mat("cell_lines/50k_TCGAbg/matrix/cells_50k.consensus.Binarymat.rds")
tcga <- load_mat("TCGA/cell_50k/matrix/TCGA_50k_cell.consensus.Binarymat.rds")

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

save(pca_res, samples, prop_var, file = "cell_tcga_pca.RData")