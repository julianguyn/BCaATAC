# load libraries
suppressPackageStartupMessages({
    library(limma)
    library(DMRcate)
})

###########################################################
# Load in data
###########################################################

# get beta values
betas <- readRDS("data/procdata/TCGA/TCGA_betas.RDS")
rownames(betas) <- betas$CpG
betas$CpG <- NULL

# read in meta data file
meta <- read.csv("metadata/TCGA_mutation_meta.csv")
meta <- meta[match(colnames(betas), meta$Sample.Name), ]

###########################################################
# Convert beta values to M values
###########################################################

# add 1e-6 to avoid log 0
mvals <- log2((betas + 1e-6) / (1 - betas + 1e-6))

###########################################################
# Differentially methylated sites
###########################################################


###########################################################
# Differentially methylated regions
###########################################################