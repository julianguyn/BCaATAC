# load libraries
suppressPackageStartupMessages({
    #library(PharmacoGx)
    #library(reshape2)
    #library(dplyr)
    #library(ComplexHeatmap)
    #library(circlize)
    library(genefu)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    #library(plyr)
    #library(ggplot2)
    #library(grid)
    #library(gridExtra)
    #library(data.table)
})

source("utils/mappings.R")
source("utils/score_bca_subtype.R")
#source("utils/palettes.R")
#source("utils/get_data.R")

###########################################################
# Load in data
###########################################################

# get PDX gene counts matrix
rna <- read.csv("data/rawdata/pdx/gene_tpm_normalized_matrix.csv")

# read in metadata (only used for checking sampleIDs)
pdx_meta <- read.csv("metadata/lupien_metadata.csv")
pdx_meta <- pdx_meta[pdx_meta$type == "PDX",]

# load in genefu subtyping models
data(pam50.robust)
data(scmgene.robust)
data(scmod2.robust)

###########################################################
# Formating RNA matrix
###########################################################

# formating sampleIDs
colnames(rna) <- sub("^S", "", colnames(rna))
colnames(rna) <- map_pdx(colnames(rna))

# add additional column for NOTCH01 passages
rna$NOTCH01P6 <- rna$NOTCH01
colnames(rna)[colnames(rna) == "NOTCH01"] <- "NOTCH01P4"
colnames(rna)[colnames(rna) == "NOTCHB01_P3_TXLCTRL1"] <- "NOTCHB01"

# check_sample_overlap(colnames(rna), pdx_meta$sampleid, "RNA", "meta")
# Samples without RNA:
# "21046" "67725" "70420" "73720"

# formating genes
rownames(rna) <- rna$X
rna$X <- NULL

###########################################################
# Map genes to create metadata
###########################################################

meta <- data.frame(
    Ensembl = rownames(rna)
)

# get gene symbol
meta$Gene.Symbol = mapIds(
    org.Hs.eg.db,
    keys=test$Ensembl,
    column="SYMBOL",
    keytype="ENSEMBL",
    multiVals="first"
)

# get entrez ids
meta$EntrezGene.ID = mapIds(
    org.Hs.eg.db,
    keys=meta$Ensembl,
    column="ENTREZID",
    keytype="ENSEMBL",
    multiVals="first"
)

###########################################################
# Compute subtyping model scores
###########################################################

pdx_pam50 <- score_bca_subtype(rna, meta, model = "pam50")
pdx_scmgene <- score_bca_subtype(rna, meta, model = "scmgene")
pdx_scmod2 <- score_bca_subtype(rna, meta, model = "scmod2")

save(
    pdx_pam50, pdx_scmgene, pdx_scmod2,
    file = "data/results/data/3-DataExploration/pdxs_subtyping_scores.RData"
)
