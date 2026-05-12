# umap - run on H4H

library(data.table)
library(ggplot2)

set.seed(101)
source("utils/palettes.R")

INDIR <- "/cluster/projects/bhklab/projects/BCaATAC/BCa_ARCHE_Scoring/data/results/"
OUTDIR <- "data/results/figures/Misc/scoring/"

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

  pca_res <- pca_res$x[,c(1:3)]
  print(head(pca_res))

  save(pca_res, samples, prop_var, file = outfile)

}

run_pca(merged_cells, "cell_tcga_pca.RData")
run_pca(merged_pdxs, "pdx_tcga_pca.RData")

###########################################################
# Load in metadata
###########################################################

# read in sample metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# read in tcga sample matching
mapping <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
mapping$Sample.Name <- gsub("\\.", "-", mapping$Sample.Name)
mapping$ATAC.Seq.File.Name <- gsub("\\.", "-", mapping$ATAC.Seq.File.Name)

###########################################################
# Helper function to plot PCA plots
###########################################################

plot_pca <- function(label) {

  # plot eigenvalues
  filename <- paste0(OUTDIR, label, "_tcga_eigen.png")
  png(filename, width = 6, height = 4, res = 600, units = "in")
  plot(prop_var,
      xlab = "Principal Components",
      ylab = "Proportion of Variance Explained",
      ylim = c(0, 1))
  dev.off()

  # add labels to pca
  pca_res <- as.data.frame(pca_res)
  for (i in 1:nrow(pca_res)) {
    if (rownames(pca_res)[i] %in% meta$filename) {
      pca_res$Sample_Type[i] <- label
      pca_res$Subtype[i] <- meta$subtype[meta$filename == rownames(pca_res)[i]]
      pca_res$Tech[i] <- meta$tech[meta$filename == rownames(pca_res)[i]]
    } else {
      pca_res$Sample_Type[i] <- "tumour"
      pca_res$Subtype[i] <- mapping$Subtype[mapping$ATAC.Seq.File.Name == rownames(pca_res)[i]]
      pca_res$Tech[i] <- "TCGA"
    }
  }
  pca_res$Subtype[is.na(pca_res$Subtype)] <- "unknown"

  # plot PC1vs2
  p <- ggplot(pca_res, aes(x = PC1, y = PC2, shape = Sample_Type, fill = Subtype)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_fill_manual(values = subtype_pal) +
    scale_shape_manual(values = c(24, 21)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  filename <- paste0(OUTDIR, label, "_tcga_pca12.png")
  ggsave(filename, p, width = 7, height = 5)

  # plot PC2vs3
  p <- ggplot(pca_res, aes(x = PC3, y = PC2, shape = Sample_Type, fill = Subtype)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_fill_manual(values = subtype_pal) +
    scale_shape_manual(values = c(24, 21)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  filename <- paste0(OUTDIR, label, "_tcga_pca23.png")
  ggsave(filename, p, width = 7, height = 5)

}

###########################################################
# Plot PCA
###########################################################

load("data/procdata/scoring/cell_tcga_pca.RData")
plot_pca("cell")

load("data/procdata/scoring/pdx_tcga_pca.RData")
plot_pca("pdx")