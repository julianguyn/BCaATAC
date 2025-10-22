setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(maftools)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  library(ggh4x)
})

source("source/MolecularSigAnalysis/helper.R")
source("source/MolecularSigAnalysis/plots.R")
source("source/palettes.R")

###########################################################
# Load in data
###########################################################

# read in meta data file
meta <- get_meta_mut()

# load in matrix file from NMF
mat <- get_ARCHE()
mat <- mat[mat$variable %in% meta$ATAC.Seq.File.Name,]

# get mafs
mafs <- merge_mafs(meta$SNV.File.Name)

###########################################################
# MAF summaries all tumours
###########################################################

plot_mafSummary(arche = "all")
plot_mafSummary(arche = "ARCHE1")
plot_mafSummary(arche = "ARCHE2")
plot_mafSummary(arche = "ARCHE3")
plot_mafSummary(arche = "ARCHE4")
plot_mafSummary(arche = "ARCHE5")
plot_mafSummary(arche = "ARCHE6")

###########################################################
# Create mutations counts matrix
###########################################################

# create counts matrix
mut <- mutCountMatrix(
  mafs,
  includeSyn = FALSE,
  countOnly = NULL,
  removeNonMutated = TRUE
)

write.table(mut, 
            file = "MolecularSigAnalysis/results/data/TCGA_mutation_matrix.tsv", 
            quote = F, sep = "\t", col.names = T, row.names = T)

###########################################################
# Map files to metadata
###########################################################

mapping <- data.frame(snv = colnames(mut))
mapping$Sample.Name <- sub("^([^-]+-[^-]+-[^-]+).*", "\\1", mapping$snv)
mapping$Sample.Name <- gsub("-", "\\.", mapping$Sample.Name)

meta$snv_label <- mapping$snv[match(meta$Sample.Name, mapping$Sample.Name)]
write.csv(meta, file = "MolecularSigAnalysis/data/TCGA_sourcefiles_mut.csv", quote = F, row.names = F)

###########################################################
# Formating dataframe for plotting
###########################################################

# filter for BCa-relevant mutations of interest
toPlot <- mut[rownames(mut) %in% bca_mutations,]

# plot presence of BCa-relevant mutations
plot_BCa_mutations(toPlot)