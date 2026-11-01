# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(patchwork)
    library(ggnewscale)
    library(ggpattern)
    library(matrixStats)
})

source("utils/get_data.R")
source("utils/palettes.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/ccl_benchmarks.R")
source("utils/plots/ARCHE_scores_heatmap.R")

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

# get normzs_cells_sumdev
zscore_cells <- get_arche_scores("data/rawdata/all_scoring/cell_tcga.Zscore.txt", c_meta)
zscore_cells_sumdev <- get_arche_sumdevs(zscore_cells, "zscore_cells", plot = TRUE)
normzs_cells_sumdev <- znorm(zscore_cells_sumdev)

# load in RNA
ubr1 <- get_pset_rna("UBR1")
ubr2 <- get_pset_rna("UBR2")
gray <- get_pset_rna("GRAY")
gcsi <- get_pset_rna("gCSI")
gdsc <- get_pset_rna("GDSC2")
ccle <- get_pset_rna("CCLE")

# load in mutation data
gdsc_mut <- get_pset_mut("GDSC2")
ccle_mut <- get_pset_mut("CCLE")

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

tozasertib_genes <- c(
    "ENSG00000087586" = "AURKA",
    "ENSG00000178999" = "AURKB",
    "ENSG00000105146" = "AURKC"
)

paclitaxel_genes <- c(
    "ENSG00000258947" = "TUBB3",
    "ENSG00000117632" = "STMN1",
    "ENSG00000047849" = "MAP4",
    "ENSG00000085563" = "ABCB1"
)

trastuzumab_genes <- c(
    "ENSG00000141736" = "ERBB2"
)

erlotinib_genes <- c(
    "ENSG00000146648" = "EGFR",
    "ENSG00000123374" = "CDK2"
)

topotecan_genes <- c(
    "ENSG00000172716" = "SLFN11"
)

PD_0325901_genes <- c(
    "ENSG00000157764" = "BRAF",
    "ENSG00000213281" = "NRAS",
    "ENSG00000133703" = "KRAS",
    "ENSG00000196712" = "NF1",
    "ENSG00000169032" = "MAP2K1"
)

# drug targets (targeted therapies are compared against target gene expression instead of PAM50)
tozasertib_target <- c("ENSG00000087586" = "AURKA")
topotecan_target <- c("ENSG00000198900" = "TOP1")
trastuzumab_target <- c("ENSG00000141736" = "ERBB2")
erlotinib_target <- c("ENSG00000146648" = "EGFR")
mek_target <- c("ENSG00000169032" = "MAP2K1") # PD-0325901 and trametinib

# ARCHE2
a2_tozasertib <- plot_rna_associations("ARCHE2", "Basal", "Tozasertib", normzs_cells_sumdev, tozasertib_genes, tozasertib_target)
a2_topotecan <- plot_rna_associations("ARCHE2", "Basal", "Topotecan", normzs_cells_sumdev, topotecan_genes, topotecan_target)

# ARCHE3
a3_trastuzumab <- plot_rna_associations("ARCHE3", "Her2", "Trastuzumab", normzs_cells_sumdev, trastuzumab_genes, trastuzumab_target)

# ARCHE4
a4_erlotinib <- plot_rna_associations("ARCHE4", "Basal", "Erlotinib", normzs_cells_sumdev, erlotinib_genes, erlotinib_target)
a4_topotecan <- plot_rna_associations("ARCHE4", "Basal", "Topotecan", normzs_cells_sumdev, topotecan_genes, topotecan_target)
a4_PD_0325901 <- plot_rna_associations("ARCHE4", "LumA", "PD-0325901", normzs_cells_sumdev, PD_0325901_genes, mek_target)
a4_trametinib <- plot_rna_associations("ARCHE4", "Basal", "Trametinib", normzs_cells_sumdev, PD_0325901_genes, mek_target)

# ARCHE5
a5_paclitaxel <- plot_rna_associations("ARCHE5", "Basal", "Paclitaxel", normzs_cells_sumdev, paclitaxel_genes) # chemotherapy: PAM50
a5_erlotinib <- plot_rna_associations("ARCHE5", "Basal", "Erlotinib", normzs_cells_sumdev, erlotinib_genes, erlotinib_target)

###########################################################
# Mutations
###########################################################

# these won't work anymore
#plot_mut_associations("ARCHE3", "Alpelisib", "PIK3CA", zscore_cells_sumdev)
#plot_mut_associations("ARCHE4", "Erlotinib", "EGFR", zscore_cells_sumdev)

genes <- c("NF1" = "ENSG00000196712", "BRAF" = "ENSG00000157764")

plot_mut_associations("ARCHE4", "Selumetinib", genes, normzs_cells_sumdev, "GDSC2", "CCLE", "CTRP")
plot_mut_associations("ARCHE4", "PD-0325901", genes, normzs_cells_sumdev, "gCSI", "GDSC2", "CCLE")
plot_mut_associations("ARCHE4", "Trametinib", genes, normzs_cells_sumdev, "GRAY", "GDSC2", "CTRP")
