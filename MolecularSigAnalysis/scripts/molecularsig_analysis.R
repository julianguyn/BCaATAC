setwd("C:/Users/julia/Documents/BCaATAC")


# load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(GSEABase)
  library(ggplot2)
  library(RColorBrewer)
  library(NMF)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(pheatmap)
  library(dplyr)
})

source("source/MolecularSigAnalysis/helper.R")
source("source/MolecularSigAnalysis/plots.R")
source("source/palettes.R")

###########################################################
# Load in data
###########################################################

# read in meta data file
meta <- read.csv("MolecularSigAnalysis/data/TCGA_sourcefiles_mut.csv")

# load in matrix file from NMF
mat <- get_ARCHE()
mat <- mat[mat$variable %in% meta$ATAC.Seq.File.Name,]

# get mafs
mafs <- merge_mafs(meta$SNV.File.Name)

# load in gene counts matrix
counts <- fread("Signatures/data/TCGA_BRCA_gene_counts.matrix") |>
  as.data.frame()
rownames(counts) <- counts$V1
counts$V1 <- NULL
counts <- t(counts)

###########################################################
# Extract BCa mutation calls
###########################################################

mafs.tnm <- trinucleotideMatrix(maf = mafs, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
mut_calls <- t(mafs.tnm$nmf_matrix)
colnames(mut_calls) = meta$Sample.Name[match(colnames(mut_calls), meta$snv_label)]
res <- list(signatures = mut_calls)

###########################################################
# Cosine Similarity against Mutational Signatures
###########################################################
 
# original 30 mutational signatures
og30.cosm <- compareSignatures(nmfRes = res, sig_db = "legacy")$cosine_similarities |> suppressMessages()
write.csv(og30.cosm, file = "MolecularSigAnalysis/results/data/og30_cosm.csv", quote = F, row.names = F)

# updated 60 signatures
v3.cosm <- compareSignatures(nmfRes = res, sig_db = "SBS")$cosine_similarities |> suppressMessages()
write.csv(v3.cosm, file = "MolecularSigAnalysis/results/data/v3_cosm.csv", quote = F, row.names = F)

###########################################################
# Single sample gsea on hallmarks gene sets
###########################################################

hm.es <- run_ssgsea(counts)

###########################################################
# Heatmaps cluster by mutational signatures
###########################################################

plot_pheatmap(og30.cosm, "Original_30")
plot_pheatmap(v3.cosm, "Updated_60")

plot_heatmap_mutsig(og30.cosm, "og30")
plot_heatmap_mutsig(v3.cosm, "v3")

###########################################################
# Correlate each pair of signatures and plot
###########################################################

# format cosine similarity matrices
og30.cosm <- t(og30.cosm) |> as.data.frame()
v3.cosm <- t(v3.cosm) |> as.data.frame()

# correlate signatures
corr_og <- corr_signatures(og30.cosm, "og30")
corr_v3 <- corr_signatures(v3.cosm, "v3")
corr_hm <- corr_signatures(hm.es, "hm")

# plot heatmaps of signature correlations
plot_molecularsig_corr(corr_og, "Original 30 COSMIC Signatures", "og30")
plot_molecularsig_corr(corr_v3, "Mutational Signatures", "v3")
plot_molecularsig_corr(corr_hm, "Hallmark Gene Sets", "hm")

# plot box plots
plot_corr_boxplots(corr_v3, corr_hm)
