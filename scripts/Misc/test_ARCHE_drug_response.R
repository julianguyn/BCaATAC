# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(ggplot2)
    library(ggh4x)
    library(reshape2)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(dplyr)
    library(readxl)
    library(data.table)
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/plots/drug_response_ccls.R")
source("utils/palettes.R")
source("utils/bca_drugs.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/plots/drug_response_pdx.R")

###########################################################
# Prepare metadata
###########################################################

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]

# remove second dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta_t <- meta[!(meta$sampleid %in% dups & meta$tech == "tina"), ]
meta_k <- meta[!(meta$sampleid %in% dups & meta$tech == "komal"), ]

###########################################################
# Load in cell line data
###########################################################

c_meta_t <- meta_t[meta_t$type == "cell_line", ]
c_meta_k <- meta_k[meta_k$type == "cell_line", ]

# load in arche zscores and normalize
zscore_50k_t <- get_arche_scores("cells", "k50", c_meta_t)
normzs_50k_t <- znorm(zscore_50k_t)
zscore_50k_k <- get_arche_scores("cells", "k50", c_meta_k)
normzs_50k_k <- znorm(zscore_50k_k)

# load in arche deviations and normalize
deviat_50k_t <- get_arche_devs("cells", "50k", c_meta_t)
normdv_50k_t <- znorm(deviat_50k_t)
deviat_50k_k <- get_arche_devs("cells", "50k", c_meta_k)
normdv_50k_k <- znorm(deviat_50k_k)

# get drug sensitivity data
load("data/procdata/CCLs/sensitivity_data.RData")

###########################################################
# Compute PC of ARCHE-drug associations in cell lines
###########################################################

# helper function to compute arche associations across all psets
arche_pc <- function(scores) {
    ubr1_PC <- compute_pc(scores, ubr1_sen, "UBR1")
    ubr2_PC <- compute_pc(scores, ubr2_sen, "UBR2")
    gray_PC <- compute_pc(scores, gray_sen, "GRAY")
    gcsi_PC <- compute_pc(scores, gcsi_sen, "gCSI")
    gdsc_PC <- compute_pc(scores, gdsc_sen, "GDSC2")
    ctrp_PC <- compute_pc(scores, ctrp_sen, "CTRP")
    ccle_PC <- compute_pc(scores, ccle_sen, "CCLE")

    # compile results
    PC_res <- rbind(ubr1_PC, ubr2_PC, gray_PC, gcsi_PC, gdsc_PC, ctrp_PC, ccle_PC)
    return(PC_res)
}

# zscores
pc_zscore_50k_t <- arche_pc(zscore_50k_t)
pc_normzs_50k_t <- arche_pc(normzs_50k_t)
pc_zscore_50k_k <- arche_pc(zscore_50k_k)
pc_normzs_50k_k <- arche_pc(normzs_50k_k)

# deviations
pc_deviat_50k_t <- arche_pc(deviat_50k_t)
pc_normdv_50k_t <- arche_pc(normdv_50k_t)
pc_deviat_50k_k <- arche_pc(deviat_50k_k)
pc_normdv_50k_k <- arche_pc(normdv_50k_k)

###########################################################
# Compile cell line results
###########################################################

pc_zscore_50k_t$Label <- "zscore_T"
pc_normzs_50k_t$Label <- "normzscr_T"
pc_zscore_50k_k$Label <- "zscore_K"
pc_normzs_50k_k$Label <- "normzscr_K"

pc_deviat_50k_t$Label <- "deviat_T"
pc_normdv_50k_t$Label <- "normdev_T"
pc_deviat_50k_k$Label <- "deviat_K"
pc_normdv_50k_k$Label <- "normdev_K"

cell_toPlot <- rbind(
    pc_zscore_50k_t, pc_normzs_50k_t,
    pc_zscore_50k_k, pc_normzs_50k_k,
    pc_deviat_50k_t, pc_normdv_50k_t,
    pc_deviat_50k_k, pc_normdv_50k_k
)
cell_toPlot$drug[cell_toPlot$drug == "945"] <- "CFI-400945"

###########################################################
# Load in PDX data
###########################################################

p_meta_t <- meta_t[meta_t$type == "PDX", ]
p_meta_k <- meta_k[meta_k$type == "PDX", ]

# load in arche zscores and normalize
zscore_t <- get_arche_scores("pdxs", "k50", p_meta_t)
normzs_t <- znorm(zscore_t)
zscore_k <- get_arche_scores("pdxs", "k50", p_meta_k)
normzs_k <- znorm(zscore_k)

# load in arche deviations and normalize
deviat_t <- get_arche_devs("PDXs", "50k", p_meta_t)
normdv_t <- znorm(deviat_t)
deviat_k <- get_arche_devs("PDXs", "50k", p_meta_k)
normdv_k <- znorm(deviat_k)

# pdx drug response data
xeva <- get_xeva("full")

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

deviat_t <- format_ARCHE(xeva, deviat_t)
normdv_t <- format_ARCHE(xeva, normdv_t)
deviat_k <- format_ARCHE(xeva, deviat_k)
normdv_k <- format_ARCHE(xeva, normdv_k)

###########################################################
# Compute PC of ARCHE-drug associations in PDXs
###########################################################

dir <- "data/results/figures/Misc/ARCHE_dr_test/"
x_zscore_t <- assess_ARCHE_PDX(zscore_t, "zscore_T", dir, plot = FALSE) |> suppressWarnings()
x_normzs_t <- assess_ARCHE_PDX(normzs_t, "normzscr_T", dir, plot = FALSE) |> suppressWarnings()
x_zscore_k <- assess_ARCHE_PDX(zscore_k, "zscore_K", dir, plot = FALSE) |> suppressWarnings()
x_normzs_k <- assess_ARCHE_PDX(normzs_k, "normzscr_K", dir, plot = FALSE) |> suppressWarnings()

x_deviat_t <- assess_ARCHE_PDX(deviat_t, "deviat_T", dir, plot = FALSE) |> suppressWarnings()
x_normdv_t <- assess_ARCHE_PDX(normdv_t, "normdev_T", dir, plot = FALSE) |> suppressWarnings()
x_deviat_k <- assess_ARCHE_PDX(deviat_k, "deviat_K", dir, plot = FALSE) |> suppressWarnings()
x_normdv_k <- assess_ARCHE_PDX(normdv_k, "normdev_K", dir, plot = FALSE) |> suppressWarnings()

###########################################################
# Compile PDX results
###########################################################

pdx_toPlot <- rbind(
    x_zscore_t, x_normzs_t,
    x_zscore_k, x_normzs_k,
    x_deviat_t, x_normdv_t,
    x_deviat_k, x_normdv_k
)
colnames(pdx_toPlot)[colnames(pdx_toPlot) == "ARCHE_label"] <- "Label"
colnames(pdx_toPlot)[colnames(pdx_toPlot) == "pairs"] <- "pair"

save(cell_toPlot, pdx_toPlot, file = "data/results/data/Misc/ARCHE_drug_response_testing.RData")

###########################################################
# Plot results
###########################################################

pair <- "ARCHE5_Paclitaxel"

# plot cells
subset <- cell_toPlot[cell_toPlot$pairs == pair,]
subset$sig <- ifelse(
    subset$FDR < 0.1,
    ifelse(subset$FDR < 0.05, 'FDR < 0.05', 'FDR < 0.1'),
    'FDR >= 0.1')
subset$text <- ifelse(
    subset$FDR < 0.1, subset$pc, NA
)

subset$Label <- factor(subset$Label, levels = c(
    "zscore_T", "zscore_K", "deviat_T", "deviat_K",
    "normzscr_T", "normzscr_K", "normdev_T", "normdev_K"
    ))

p <- ggplot(subset, aes(x = pset, y = Label, fill = pc, shape = sig, size = -log(FDR))) +
    geom_point() +
    geom_text(aes(label = round(text, 2)), color = "black", size = 2.5) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_size(range = c(2, 12)) +
    scale_fill_gradient2(
        low = "#BC4749",
        high = "#689CB0",
        mid = "#C2BBC9",
        limits = c(-0.8, 0.8)
    ) +
    theme_bw() +
    theme(legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size = 10)) +
    ggtitle(pair)

filename <- paste0("data/results/figures/Misc/ARCHE_dr_test/cells_", pair, ".png")
ggsave(filename, p, w = 6, h = 4)

