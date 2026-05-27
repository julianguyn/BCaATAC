# PCA - run on H4H
# cd /cluster/projects/bhklab/projects/BCaATAC/BCa_ARCHE_Scoring/temp_scripts

library(data.table)

set.seed(101)


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
cells <- load_mat("cell_lines/50k_TCGAbg/matrix/cells_50k.consensus.Binarymat.rds")
tcga_cells <- load_mat("TCGA/cell_50k/matrix/TCGA_50k_cell.consensus.Binarymat.rds")

# load in pdxs
pdxs <- load_mat("PDXs/50k_TCGAbg/matrix/PDXs_50k.consensus.Binarymat.rds")
tcga_pdxs <- load_mat("TCGA/pdx_50k/matrix/TCGA_50k_pdx.consensus.Binarymat.rds")

###########################################################
# Merge matrices and format for pca
###########################################################

# helper function to format dataframes for pca
format_df <- function(sample, tcga) {

  # merged on peaks
  merged <- merge(
    sample,
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

saveRDS(merged_cells, file = "merged_cells.rds") #1923648 peaks
saveRDS(merged_pdxs, file = "merged_pdxs.rds") #1869335 peaks

###########################################################
# Identify differential peaks
###########################################################
# -- run on H4H

# load in RDSs
print("loading in data")
merged_cells <- readRDS("merged_cells.rds")
merged_pdxs <- readRDS("merged_pdxs.rds")

find_diff_peaks <- function(merged, n_samples) {
  
  samples <- merged[1:n_samples,]
  tumours <- merged[(n_samples+1):nrow(merged),]

  sample_count <- colSums(samples)
  sample_zero <- sample_count[sample_count == 0]
  sample_ones <- sample_count[sample_count == n_samples]

  tumour_count <- colSums(tumours)
  tumour_zero <- tumour_count[tumour_count == 0]
  tumour_ones <- tumour_count[tumour_count == 75]

  diff_peaks <- c(
    intersect(names(sample_zero), names(tumour_ones)),
    intersect(names(sample_ones), names(tumour_zero))
  )
  return(diff_peaks)
}

cell_diff_peaks <- find_diff_peaks(merged_cells, 64)
saveRDS(cell_diff_peaks, file = "results/cell_diff_peaks.rds")
pdx_diff_peaks <- find_diff_peaks(merged_pdxs, 88)
saveRDS(pdx_diff_peaks, file = "results/pdx_diff_peaks.rds")

###########################################################
# Fisher's test to find differentially accessible peaks
###########################################################
# -- run on H4H


# sample groups
cell_groups <- factor(c(rep("Cells", 64), rep("Tumour", 75)))
pdx_groups <- factor(c(rep("PDXs", 88), rep("Tumour", 75)))


run_fisher <- function(merged, group) {

  merged <- as.data.frame(merged)
  # extract peak names

  results <- lapply(colnames(merged), function(peak) {
    peak_vals <- merged[[peak]]
    tbl <- table(peak_vals, group)
    
    #catch cases where peak is all 0 or 1 in one group
    if (nrow(tbl) < 2) {
      data.frame(
        peak = peak,
        ft = FALSE,
        p.value = NA,
        odds_ratio = NA,
        row.names = NULL
      )
    } else {
      ft <- fisher.test(tbl)
      data.frame(
        peak = peak,
        ft = TRUE,
        p.value = ft$p.value,
        odds_ratio = as.numeric(ft$estimate),
        row.names = NULL
      )
    }
  })

  results_df <- do.call(rbind, Filter(Negate(is.null), results))
  results_df$padj <- p.adjust(results_df$p.value, method = "BH")
  results_df <- results_df[order(results_df$padj), ]

}

print("running cell line")
cell_res <- run_fisher(merged_cells, cell_groups)
print("running PDXs")
pdx_res <- run_fisher(merged_pdxs, pdx_groups)

print("saving")
saveRDS(cell_res, file = "results/cells_fisher.rds")
saveRDS(pdx_res, file = "results/pdxs_fisher.rds")

###########################################################
# Run PCA
###########################################################

library(ggplot2)
source("utils/palettes.R")

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

  pca_res <- pca_res$x[,c(1:5)]
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

# helper function to plot PCA
helper_plot_pca <- function(pca_res, PCx, PCy, label) {

  nums <- paste0(sub("PC", "", PCx), sub("PC", "", PCy))
  
  p <- ggplot(pca_res, aes(x = .data[[PCx]], y = .data[[PCy]], shape = Sample_Type, fill = Subtype)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_fill_manual(values = subtype_pal) +
    scale_shape_manual(values = c(24, 21)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
  filename <- paste0(OUTDIR, label, "_tcga_pca", nums, ".png")
  ggsave(filename, p, width = 7, height = 5)

}

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

  # plot PCs
  helper_plot_pca(pca_res, "PC1", "PC2", label)
  helper_plot_pca(pca_res, "PC2", "PC3", label)
  helper_plot_pca(pca_res, "PC3", "PC4", label)
  helper_plot_pca(pca_res, "PC4", "PC5", label)
  helper_plot_pca(pca_res, "PC5", "PC6", label)
  helper_plot_pca(pca_res, "PC6", "PC7", label)
  helper_plot_pca(pca_res, "PC7", "PC8", label)
  helper_plot_pca(pca_res, "PC8", "PC9", label)
  helper_plot_pca(pca_res, "PC9", "PC10", label)
}

###########################################################
# Plot PCA
###########################################################

load("data/procdata/scoring/cell_tcga_pca.RData")
plot_pca("cell")

load("data/procdata/scoring/pdx_tcga_pca.RData")
plot_pca("pdx")
