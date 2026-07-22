# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(matrixStats)
    library(ggh4x)
    library(reshape2)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(dplyr)
    library(ROCR)
    library(readxl)
    library(data.table)
    library(patchwork)
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/plots/drug_response_ccls.R")
source("utils/palettes.R")
source("utils/bca_drugs.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/plots/drug_response_pdx.R")
source("utils/plots/drug_response_pdx_indivplots.R")


###########################################################
# Prepare metadata
###########################################################

# read in sample metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]

# remove second dups
dups <- meta$sampleid[duplicated(meta$sampleid)] # only PDXs have dups
meta_t <- meta[!(meta$sampleid %in% dups & meta$tech == "tina"), ]
meta_k <- meta[!(meta$sampleid %in% dups & meta$tech == "komal"), ]

###########################################################
# Load in PDX data
###########################################################

p_meta_t <- meta_t[meta_t$type == "PDX", ]
p_meta_k <- meta_k[meta_k$type == "PDX", ]

# load in arche zscores and normalize
zscore_t <- get_arche_scores("data/rawdata/all_scoring/pdx_tcga.Zscore.txt", p_meta_t)
normzs_t <- znorm(zscore_t)
zscore_k <- get_arche_scores("data/rawdata/all_scoring/pdx_tcga.Zscore.txt", p_meta_k)
normzs_k <- znorm(zscore_k)

# pdx drug response data
xeva <- get_xeva("full")

# read in RNA
rna <- read.csv("data/rawdata/pdx/gene_tpm_normalized_matrix.csv")
rownames(rna) <- rna$X
colnames(rna) <- sub("^S", "", colnames(rna))
colnames(rna) <- map_pdx(colnames(rna))
rna <- rna[,colnames(rna) %in% p_meta_t$sampleid]

###########################################################
# Removing low deviation samples
###########################################################

# get samples above threshold
zscore_t_sumdev <- get_arche_sumdevs(zscore_t, "zscore_pdx_t", plot = TRUE)
zscore_k_sumdev <- get_arche_sumdevs(zscore_k, "zscore_pdx_k", plot = TRUE)

# normalize
normzs_t_sumdev <- znorm(zscore_t_sumdev)
normzs_k_sumdev <- znorm(zscore_k_sumdev)

###########################################################
# Assign ARCHE scores to PDX xeva
###########################################################

# helper function
format_ARCHE <- function(xeva, scores) {
    scores <- t(scores) |> as.data.frame()
    df <- xeva[xeva$patient.id %in% rownames(scores),]
    df$ARCHE6 <- df$ARCHE5 <- df$ARCHE4 <- df$ARCHE3 <- df$ARCHE2 <- df$ARCHE1 <- NA
    for (i in 1:nrow(df)) {
        sample <- df$patient.id[i]
        df$ARCHE1[i] <- scores[rownames(scores) == sample,]$ARCHE1
        df$ARCHE2[i] <- scores[rownames(scores) == sample,]$ARCHE2
        df$ARCHE3[i] <- scores[rownames(scores) == sample,]$ARCHE3
        df$ARCHE4[i] <- scores[rownames(scores) == sample,]$ARCHE4
        df$ARCHE5[i] <- scores[rownames(scores) == sample,]$ARCHE5
        df$ARCHE6[i] <- scores[rownames(scores) == sample,]$ARCHE6
        
    }
    return(df)
}

# get ARCHE scores
zscore_t <- format_ARCHE(xeva, zscore_t)
normzs_t <- format_ARCHE(xeva, normzs_t)
zscore_k <- format_ARCHE(xeva, zscore_k)
normzs_k <- format_ARCHE(xeva, normzs_k)

zscore_t_sumdev <- format_ARCHE(xeva, zscore_t_sumdev)
normzs_t_sumdev <- format_ARCHE(xeva, normzs_t_sumdev)
zscore_k_sumdev <- format_ARCHE(xeva, zscore_k_sumdev)
normzs_k_sumdev <- format_ARCHE(xeva, normzs_k_sumdev)

###########################################################
# Compute PC of ARCHE-drug associations in PDXs
###########################################################

plot_dir <- "indiv"
x_zscore_t <- assess_ARCHE_PDX(zscore_t, "zscore_T", plot_dir, plot = FALSE) |> suppressWarnings()
x_normzs_t <- assess_ARCHE_PDX(normzs_t, "normzscr_T", plot_dir, plot = FALSE) |> suppressWarnings()
x_zscore_k <- assess_ARCHE_PDX(zscore_k, "zscore_K", plot_dir, plot = FALSE) |> suppressWarnings()
x_normzs_k <- assess_ARCHE_PDX(normzs_k, "normzscr_K", plot_dir, plot = FALSE) |> suppressWarnings()

x_zscore_t_sumdev <- assess_ARCHE_PDX(zscore_t_sumdev, "zscore_sumdev_T", plot_dir, plot = FALSE) |> suppressWarnings()
x_normzs_t_sumdev <- assess_ARCHE_PDX(normzs_t_sumdev, "normzscr_sumdev_T", plot_dir, plot = FALSE) |> suppressWarnings()
x_zscore_k_sumdev <- assess_ARCHE_PDX(zscore_k_sumdev, "zscore_sumdev_K", plot_dir, plot = FALSE) |> suppressWarnings()
x_normzs_k_sumdev <- assess_ARCHE_PDX(normzs_k_sumdev, "normzscr_sumdev_K", plot_dir, plot = FALSE) |> suppressWarnings()

###########################################################
# Compile PDX results
###########################################################

# all samples
pdx_all_toPlot <- rbind(
    x_zscore_t, x_normzs_t,
    x_zscore_k, x_normzs_k,
    x_zscore_t_sumdev, x_normzs_t_sumdev,
    x_zscore_k_sumdev, x_normzs_k_sumdev
)
colnames(pdx_all_toPlot)[colnames(pdx_all_toPlot) == "ARCHE_label"] <- "Label"
colnames(pdx_all_toPlot)[colnames(pdx_all_toPlot) == "pair"] <- "pairs"

save(pdx_all_toPlot,
    file = "data/results/data/4-DrugResponse/PDX/pdx_response_july222026.RData")

###########################################################
# Plot compiled PDX results
###########################################################

plot_compiled_PDX(pdx_all_toPlot)

###########################################################
# Assess ADCs
###########################################################

plot_all_ADC_ARCHE <- function(drug, gene) {

    plot_ADC_ARCHE(zscore_t, "zscore_T", drug, gene)
    plot_ADC_ARCHE(zscore_k, "zscore_K", drug, gene)
    plot_ADC_ARCHE(normzs_t, "normzscr_T", drug, gene)
    plot_ADC_ARCHE(normzs_k, "normzscr_K", drug, gene)

    plot_ADC_ARCHE(zscore_t_sumdev, "zscore_sumdev_T", drug, gene)
    plot_ADC_ARCHE(zscore_k_sumdev, "zscore_sumdev_K", drug, gene)
    plot_ADC_ARCHE(normzs_t_sumdev, "normzscr_sumdev_T", drug, gene)
    plot_ADC_ARCHE(normzs_k_sumdev, "normzscr_sumdev_K", drug, gene)
}

plot_all_ADC_ARCHE("SACITUZAMAB-GOVITECAN", "ENSG00000184292")
plot_all_ADC_ARCHE("TRASTUZUMAB-DERUXTECAN", "ENSG00000141736")
plot_all_ADC_ARCHE("DATOPOTAMAB-DERUXTECAN", "ENSG00000184292")
plot_all_ADC_ARCHE("AZD-8205", "ENSG00000134258")
