# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
})

source("utils/get_data.R")
source("utils/palettes.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/ccl_benchmarks.R")

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

pam50_scores <- rbind(
    format_pam50(ubr1_pam50, "UBR1"), format_pam50(ubr2_pam50, "UBR2"), format_pam50(gray_pam50, "GRAY"),
    format_pam50(gcsi_pam50, "gCSI"), format_pam50(gdsc_pam50, "GDSC2"), format_pam50(ccle_pam50, "CCLE"),
    format_pam50(ccle_pam50, "CTRP")
)

###########################################################
# Plot associations
###########################################################

# compile gene expression
paclitaxel_genes <- c(
    "ENSG00000258947" = "TUBB3",
    "ENSG00000117632" = "STMN1",
    "ENSG00000047849" = "MAP4",
    "ENSG00000085563" = "ABCB1"
)

trastuzumab_genes <- c(
    "ENSG00000141736" = "ERBB2"
)

# ARCHE3
a3_trastuzumab <- plot_associations("ARCHE3", "Her2", "Trastuzumab", zscore_cells_sumdev, trastuzumab_genes)

# ARCHE5
a5_paclitaxel <- plot_associations("ARCHE5", "Basal", "Paclitaxel", zscore_cells_sumdev, paclitaxel_genes)
