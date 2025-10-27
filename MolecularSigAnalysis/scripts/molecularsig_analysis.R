# load libraries
suppressPackageStartupMessages({
  library(maftools)
  library(data.table)
  library(GSVA)
  library(GSEABase)
  library(ggplot2)
  library(RColorBrewer)
  library(NMF)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(pheatmap)
  library(dplyr)
  library(qusage)
})

source("source/MolecularSigAnalysis/helper.R")
source("source/MolecularSigAnalysis/plots.R")
source("source/DataExploration/helper.R")
source("source/DrugResponse/helper.R")
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

# load in tumour gene counts matrix
t_counts <- fread("Signatures/data/TCGA_BRCA_gene_counts.matrix") |>
  as.data.frame()
rownames(t_counts) <- t_counts$V1
t_counts$V1 <- NULL
t_counts <- t(t_counts)

# load in ccls gene counts matrix
c_counts <- get_ubr2_rna()
colnames(c_counts) <- gsub("\\.", "-", colnames(c_counts))

# quick cell line mapping from map_sen()
for (i in 1:ncol(c_counts)) {
    cell = colnames(c_counts)[i]
    if (cell %in% names(mapping_cells)) {colnames(c_counts)[i] <- unname(mapping_cells[cell])}
}

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

# load in gmt files
hallmarks <- getGmt("MolecularSigAnalysis/data/h.all.v2025.1.Hs.symbols.gmt")  # from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
myc_targs <- read.gmt("MolecularSigAnalysis/data/All_MYC_Target_Signatures.gmt") # from peter lin

# ssgsea on tumour counts
t_hm.es <- run_ssgsea(t_counts, hallmarks, "hallmarks")
t_my.es <- run_ssgsea(t_counts, myc_targs, "MYC")

# ssgsea on ccls counts
c_hm.es <- run_ssgsea(c_counts, hallmarks, "hallmarks_ccls")
c_my.es <- run_ssgsea(c_counts, myc_targs, "MYC_ccls")

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
corr_t_hm <- corr_signatures(t_hm.es, "hm")
corr_t_my <- corr_signatures(t_my.es, "myc")
corr_c_hm <- corr_signatures(c_hm.es, "hm_ccls", "ccls")
corr_c_my <- corr_signatures(c_my.es, "myc_ccls", "ccls")

# plot heatmaps of signature correlations (tumours)
plot_molecularsig_corr(corr_og, "Original 30 COSMIC Signatures", "og30")
plot_molecularsig_corr(corr_v3, "Mutational Signatures", "v3")
plot_molecularsig_corr(corr_t_hm, "Hallmark Gene Sets", "hm")
plot_molecularsig_corr(corr_t_my, "MYC Target Signatures", "myc")
plot_molecularsig_corr(corr_c_hm, "Hallmark Gene Sets", "hm", "ccls")
plot_molecularsig_corr(corr_c_my, "MYC Target Signatures", "myc", "ccls")

# plot box plots
plot_corr_boxplots(corr_v3, corr_t_hm)
