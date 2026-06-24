suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
    library(tidyverse)
    library(patchwork)
    library(genefu)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
})

source("utils/get_data.R")
source("utils/palettes.R")

set.seed(123)

###########################################################
# Load in data
###########################################################

# load in RNA matrix
rna <- get_tcga_rna()
colnames(rna) <- gsub("\\.", "-", colnames(rna))
rna <- t(rna)

# load in genefu subtyping models
data(pam50.robust)
data(scmgene.robust)
data(scmod2.robust)

###########################################################
# Build gene annotation
###########################################################

gene_symbols <- colnames(rna)

entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = gene_symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
)

annot <- data.frame(
  probe = gene_symbols,
  EntrezGene.ID = as.integer(entrez_ids),
  row.names = gene_symbols,
  stringsAsFactors = FALSE
)

###########################################################
# PAM50 subtyping
###########################################################

pam50_res <- molecular.subtyping(
  sbt.model = "pam50",
  data = rna,
  annot = annot,
  do.mapping = TRUE
)

saveRDS(pam50_res, file = "data/procdata/TCGA/pam50_subtyping.rds")