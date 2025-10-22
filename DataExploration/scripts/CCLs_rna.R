# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
})

source("source/DataExploration/helper.R")
source("source/DrugResponse/helper.R")

###########################################################
# Load in data
###########################################################

# get UBR2 gene counts
ubr2 <- load_bca_RNA()