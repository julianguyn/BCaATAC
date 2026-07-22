# test ARCHE scores with TCGA peaks as the bakground

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

plot_dir <- paste0("data/results/figures/Misc/all_scoring/pdx_tcga")
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

# last save: June 18, 2026
#save(cell_toPlot, pdx_all_toPlot, pdx_sumdev_toPlot,
#    file = "data/results/data/Misc/sample_tcga_ARCHE_drug_response_testing_all_scoring.RData")


###########################################################
# Plot PDX results
###########################################################

compile <- pdx_all_toPlot

sig <- compile[which(abs(compile$PC.BAR_median) > 0.4 & compile$pval.BAR_median < 0.1),]
sig_pairs <- sig$pair
toPlot <- compile[compile$pair %in% sig_pairs,]
toPlot$sig <- ifelse(toPlot$pval.BAR_median < 0.1, 'pval < 0.1', 'pval >= 0.1')

toPlot$Label <- factor(toPlot$Label, levels = c(
    "zscore_T", "zscore_K", "normzscr_T", "normzscr_K",
    "zscore_sumdev_T", "zscore_sumdev_K", "normzscr_sumdev_T", "normzscr_sumdev_K"
    )
)

for (arche in paste0("ARCHE", 1:6)) {

    subset <- toPlot[toPlot$ARCHE == arche,]
    subset$score <- subset$Label
    #subset$score <- ifelse(
    #    sub("_.*", "", subset$Label) %in% c("zscore", "normzscr"),
    #    ifelse(sub(".*_", "", subset$Label) == "T", "ZScore_T", "ZScore_K"),
    #    ifelse(sub(".*_", "", subset$Label) == "T", "Dev_T", "Dev_K"))
    

    p1 <- ggplot(subset, aes(x = drug, y = score, fill = N)) +
        geom_tile() +
        geom_text(aes(label = N), size = 2.5) +
        scale_fill_gradient(high = "#B6B8D6", low = "#BBDBD1") +
        theme_void() +
        theme(
            axis.text.y = element_text(size = 8, hjust = 1),
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
            axis.text.x = element_text(size=6, angle=25, hjust=1, vjust=1, margin = margin(t = 3))
        ) 
    
    p <- p1 / p2 + plot_layout(heights = c(3, 6))
    filename <- paste0("data/results/figures/4-DrugResponse/PDX/", arche, ".png")
    ggsave(filename, p, w = 10, h = 5)
}
