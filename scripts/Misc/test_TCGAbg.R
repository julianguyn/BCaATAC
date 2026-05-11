# umap - run on H4H

library(data.table)
library(umap)

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

# format for umap
merged <- as.data.frame(merged)
rownames(merged) <- paste0(merged$seqnames, ":", merged$start, ":", merged$end)
merged <- merged[,-c(1:3)]
merged <- t(merged)

umap_res <- umap(merged)
umap_res <- umap_res$layout |> as.data.frame()
umap_res$sample <- rownames(merged)

save(umap_res, file = "cell_tcga_umap.RData")