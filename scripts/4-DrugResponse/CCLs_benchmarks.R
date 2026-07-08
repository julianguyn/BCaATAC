# load libraries
suppressPackageStartupMessages({
    library(data.table)


    library(PharmacoGx)
    library(survcomp)
    library(wesanderson)
    library(ggplot2)
    library(ggh4x)
    library(reshape2)
    library(meta)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(dplyr)
    library(readxl)
})

source("utils/get_data.R")
source("utils/palettes.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")

###########################################################
# Prepare metadata
###########################################################

# read in sample metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]

# get cell lines
c_meta <- meta[meta$type == "cell_line", ]

###########################################################
# Load in cell line data
###########################################################

# sumdev zscores
zscore_cells <- get_arche_scores(paste0("data/rawdata/all_scoring/cell_tcga.Zscore.txt"), c_meta)
zscore_cells_sumdev <- get_arche_sumdevs(zscore_cells)

# load in RNA
ubr1 <- get_pset_rna("UBR1")
ubr2 <- get_pset_rna("UBR2")
gray <- get_pset_rna("GRAY")
gcsi <- get_pset_rna("gCSI")
gdsc <- get_pset_rna("GDSC2")
ccle <- get_pset_rna("CCLE")

# load in PAM50 subtyping
load("data/results/data/3-DataExploration/ccls_subtyping_scores.RData")

# get drug sensitivity data
load("data/procdata/CCLs/sensitivity_data.RData")


###########################################################
# Compile PAM50 data
###########################################################

format_pam50 <- function(pam50, label) {
    pam50$Label <- label
    pam50$Sample <- rownames(pam50)
    rownames(pam50) <- NULL
    return(pam50)
}

pam50_scores <- rbind(
    format_pam50(ubr1_pam50, "UBR1"), format_pam50(ubr2_pam50, "UBR2"), format_pam50(gray_pam50, "GRAY"),
    format_pam50(gcsi_pam50, "gCSI"), format_pam50(gdsc_pam50, "GDSC2"), format_pam50(ccle_pam50, "CCLE")
)

###########################################################
# Indiv plots for associations of interest
###########################################################

# create dataframe to store results
pcc <- data.frame(matrix(nrow=0, ncol=5))
colnames(pcc) <- c("Feature_Drug", "Label", "PSet", "PCC", "pvalue")

# plot individual plots
plot_indivPlot("ARCHE5_Paclitaxel", zscore_cells_sumdev, "ARCHE", c_meta)
plot_indivPlot("Basal_Paclitaxel", pam50_scores, "PAM50", c_meta)