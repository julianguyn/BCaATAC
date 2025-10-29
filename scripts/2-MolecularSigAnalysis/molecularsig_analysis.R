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
  library(ggpubr)
  library(grid)
  library(gridExtra)
  library(ggh4x)
})

source("utils/plots/molecularsig_analysis.R")
source("utils/corr_signatures.R")
source("utils/get_data.R")
source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# read in meta data file
meta <- read.csv("metadata/TCGA_mutation_meta.csv")

# load in matrix file from NMF
mat <- get_arche_tcga()
mat <- mat[mat$variable %in% meta$ATAC.Seq.File.Name,]

# get mafs
mafs <- merge_mafs(meta$SNV.File.Name)

# load in tumour gene counts matrix
t_counts <- get_tcga_rna()

# load in ccls gene counts matrix
c_counts <- get_pset_rna("UBR2", gene.symbol = TRUE) |> as.matrix()

# load in gmt files
hallmarks <- getGmt("data/rawdata/gmt/h.all.v2025.1.Hs.symbols.gmt")  
myc_targs <- read.gmt("data/rawdata/gmt/All_MYC_Target_Signatures.gmt") 

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

# run cosine similarity
og30.cosm <- compareSignatures(nmfRes = res, sig_db = "legacy")$cosine_similarities |> suppressMessages()
v3.cosm <- compareSignatures(nmfRes = res, sig_db = "SBS")$cosine_similarities |> suppressMessages()

# save results
write.csv(og30.cosm, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/og30_cosm.csv", quote = F, row.names = F)
write.csv(v3.cosm, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/v3_cosm.csv", quote = F, row.names = F)

###########################################################
# Single sample gsea on hallmarks gene sets
###########################################################

# ssgsea on tumour counts
t_hm.es <- gsva(ssgseaParam(t_counts, hallmarks), verbose = FALSE)
t_my.es <- gsva(ssgseaParam(t_counts, myc_targs), verbose = FALSE)

# ssgsea on ccls counts
c_hm.es <- gsva(ssgseaParam(c_counts, hallmarks), verbose = FALSE)
c_my.es <- gsva(ssgseaParam(c_counts, myc_targs), verbose = FALSE)

# save results
write.csv(t_hm.es, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/hallmarks_ES.csv", quote = F, row.names = T)
write.csv(t_my.es, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/MYC_ES.csv", quote = F, row.names = T)
write.csv(c_hm.es, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/hallmarks_ccls_ES.csv", quote = F, row.names = T)
write.csv(c_my.es, file = "data/results/data/2-MolecularSigAnalysis/molecularsig/MYC_ccls_ES.csv", quote = F, row.names = T)

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
