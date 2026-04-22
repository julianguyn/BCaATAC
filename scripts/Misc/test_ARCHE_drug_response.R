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
# Try removing low deviation samples
###########################################################

# helper function to get sum of magnitude deviations
get_devs <- function(df, label, lim1, lim2) {

    devs <- colSums(abs(df)) |> as.data.frame()
    colnames(devs) <- "Sum"
    n1 <- paste0(as.character(length(devs[devs$Sum > lim1,])), "/", nrow(devs))
    n2 <- paste0(as.character(length(devs[devs$Sum > lim2,])), "/", nrow(devs))

    p <- ggplot(devs, aes(x = Sum)) +
        geom_histogram(fill = random_lightblue, color = "black", linewidth = 0.3) +
        geom_text(stat = "bin",
            aes(label = after_stat(count)),
            vjust = -0.5, size = 3) +
        theme_bw() +
        ggtitle(paste0(label, ";  n >", lim1, ":", n1, ";  n >", lim2, ":", n2))
    filename <- paste0("data/results/figures/Misc/sumdevs/", label, ".png")
    ggsave(filename, p, w=5, h=4)

    to_keep <- rownames(devs[devs$Sum > lim1,])
    return(to_keep)
}

zscore_50k_t_samples <- get_devs(zscore_50k_t, "zscore_T", 150, 200)
zscore_50k_k_samples <- get_devs(zscore_50k_k, "zscore_K", 150, 200)
deviat_50k_t_samples <- get_devs(deviat_50k_t, "deviat_T", 0.5, 0.6)
get_devs(deviat_50k_k, "deviat_K", 0.5, 0.6)

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
colnames(pdx_toPlot)[colnames(pdx_toPlot) == "pair"] <- "pairs"

save(cell_toPlot, pdx_toPlot, file = "data/results/data/Misc/ARCHE_drug_response_testing.RData")

###########################################################
# Plot PDX results
###########################################################

for (arche in paste0("ARCHE", 1:6)) {
    compile <- pdx_toPlot
    sig <- compile[which(abs(compile$PC.BAR_median) > 0.4 & compile$pval.BAR_median < 0.1),]
    sig_pairs <- sig$pair
    toPlot <- compile[compile$pair %in% sig_pairs,]
    toPlot$sig <- ifelse(toPlot$pval.BAR_median < 0.1, 'pval < 0.1', 'pval >= 0.1')

    toPlot$Label <- factor(toPlot$Label, levels = c(
        "zscore_T", "zscore_K", "deviat_T", "deviat_K",
        "normzscr_T", "normzscr_K", "normdev_T", "normdev_K")
    )

    subset <- toPlot[toPlot$ARCHE == arche,]

    p1 <- ggplot(subset, aes(x = drug, y = "No. Samples", fill = N)) +
        geom_tile() +
        geom_text(data = subset, aes(label = N)) +
        scale_fill_gradient(high = "#B6B8D6", low = "#BBDBD1") +
        theme_void() +
        theme(
            legend.position = "none") +
        ggtitle(arche)

    p2 <- ggplot(subset, aes(x = drug, y = Label, fill = PC.BAR_median, size = -log(pval.BAR_median), shape = sig)) +
        geom_point() +
        geom_text(data = subset(subset, sig == 'pval < 0.1'), aes(label = round(PC.BAR_median, 2)), color = "black", size = 2.5) +
        scale_shape_manual(values = c(21, 24)) +
        scale_size(range = c(2, 12)) +
        scale_fill_gradient2(
            low = "#BC4749",
            high = "#689CB0",
            mid = "#C2BBC9",
            limits = c(-1, 1)
        ) +
        theme_bw() +
        theme(
            legend.key.size = unit(0.3, 'cm'),
            axis.text.x = element_text(size=6, angle=45, hjust=1, vjust=1, margin = margin(t = 3))
        ) 
    
    p <- p1 / p2 + plot_layout(heights = c(1, 6))
    filename <- paste0("data/results/figures/Misc/ARCHE_dr_test/PDXs_", arche, ".png")
    ggsave(filename, p, w = 10, h = 5)
}

###########################################################
# Plot cell results
###########################################################

# plot cells
plot_cells <- function(pair) {
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
}

plot_cells("ARCHE5_Paclitaxel")