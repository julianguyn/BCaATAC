# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(data.table)
    library(tidyverse)
    library(ggpubr)
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# load subtyping scores
load("data/results/data/3-DataExploration/ccls_subtyping_scores.RData")
gray_pam50$Subtype <- NULL

# get drug sensitivity data
load("data/procdata/CCLs/sensitivity_data.RData")
gray_sen <- t(gray_sen) |> as.data.frame()

# load in RNA expression
gray_rna <- get_pset_rna("GRAY")

# function to scale to 1e6
scale_rna <- function(x) {
  x / sum(x, na.rm = TRUE) * 1e6
}
gray_rna <- apply(gray_rna, 2, scale_rna) |> as.data.frame()
gray_rna <- t(gray_rna) |> as.data.frame()

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]
c_meta <- meta[meta$type == "cell_line", ]

# load in arche scores
cells_50k <- get_arche_scores("cells", "k50", c_meta)
norm_50k <- t(znorm(cells_50k)) |> as.data.frame()
cells_50k <- t(cells_50k) |> as.data.frame()

###########################################################
# Function to plot correlations
###########################################################

plot_association <- function(df, drug, variable) {

    # get correlation
    res <- cor.test(df[[drug]], df[[variable]], method = "pearson")
    pcc <- round(res$estimate, 4)
    pval <- round(res$p.value, 5)

    # make plot
    p <- ggplot() +
        geom_smooth(
            data = df,
            aes(x = .data[[variable]], y = .data[[drug]]),
            method = "lm",
            color = "gray"
        ) +
        geom_point(
            data = df,
            aes(x = .data[[variable]], y = .data[[drug]], fill = subtype),
            shape = 21, size = 3
        ) + 
        scale_fill_manual("Assigned\nBCa Subtype", values = subtype_pal) +
        theme_minimal() +
        theme(panel.border = element_rect(), legend.key.size = unit(0.5, "cm"),) +
        ggtitle(paste("PCC:", pcc, ", pval:", pval)) +
        labs(y = paste(drug, "Response (AAC)"))
    return(p)
}

# make dataframe
df <- data.frame(Sample = unique(c(rownames(gray_pam50), rownames(gray_rna), rownames(gray_sen))))
df <- df[order(df$Sample, decreasing = FALSE),,drop=FALSE]
df$subtype <- c_meta$subtype[match(df$Sample, meta$sampleid)]
df$ARCHE1 <- cells_50k$ARCHE1[match(df$Sample, rownames(cells_50k))]
df$ARCHE2 <- cells_50k$ARCHE2[match(df$Sample, rownames(cells_50k))]
df$ARCHE3 <- cells_50k$ARCHE3[match(df$Sample, rownames(cells_50k))]
df$ARCHE4 <- cells_50k$ARCHE4[match(df$Sample, rownames(cells_50k))]
df$ARCHE5 <- cells_50k$ARCHE5[match(df$Sample, rownames(cells_50k))]
df$ARCHE6 <- cells_50k$ARCHE6[match(df$Sample, rownames(cells_50k))]

###########################################################
# Paclitaxel - all cells
###########################################################

toPlot <- df

# get variables
toPlot$Basal <- gray_pam50$Basal[match(toPlot$Sample, rownames(gray_pam50))]
toPlot$Paclitaxel <- gray_sen$Paclitaxel[match(toPlot$Sample, rownames(gray_sen))]

p1 <- plot_association(toPlot, "Paclitaxel", "Basal")
p2 <- plot_association(toPlot, "Paclitaxel", "ARCHE5")

png("data/results/figures/Misc/drugresponse/paclitaxel.png", width=7, height=2.5, units='in', res = 600, pointsize=80)
ggarrange(p1, p2, common.legend = TRUE, legend = "right")
dev.off()


###########################################################
# Gsk-461364 - all cells
###########################################################

toPlot <- df

# get variables
toPlot$Basal <- gray_pam50$Basal[match(toPlot$Sample, rownames(gray_pam50))]
toPlot$PLK1 <- gray_rna$ENSG00000166851[match(toPlot$Sample, rownames(gray_rna))]
toPlot$'Gsk-461364' <- gray_sen$'Gsk-461364'[match(toPlot$Sample, rownames(gray_sen))]

#p1 <- plot_association(toPlot, "Gsk-461364", "Basal")
p1 <- plot_association(toPlot, "Gsk-461364", "PLK1")
p2 <- plot_association(toPlot, "Gsk-461364", "ARCHE5")

png("data/results/figures/Misc/drugresponse/gsk.png", width=7, height=2.5, units='in', res = 600, pointsize=80)
ggarrange(p1, p2, common.legend = TRUE, legend = "right")
dev.off()


###########################################################
# Etoposide - all cells
###########################################################

toPlot <- df

# get variables
toPlot$LumB <- gray_pam50$LumB[match(toPlot$Sample, rownames(gray_pam50))]
toPlot$TOP2A <- gray_rna$ENSG00000131747[match(toPlot$Sample, rownames(gray_rna))]
toPlot$Etoposide <- gray_sen$Etoposide[match(toPlot$Sample, rownames(gray_sen))]

#p1 <- plot_association(toPlot, "Etoposide", "LumB")
p1 <- plot_association(toPlot, "Etoposide", "TOP2A")
p2 <- plot_association(toPlot, "Etoposide", "ARCHE4")

png("data/results/figures/Misc/drugresponse/etoposide.png", width=7, height=2.5, units='in', res = 600, pointsize=80)
ggarrange(p1, p2, common.legend = TRUE, legend = "right")
dev.off()


#####################################################################################
#####################################################################################
# Make PDX plots smaller to fit on page lol
#####################################################################################
#####################################################################################

source("utils/plots/drug_response_pdx_indivplots.R")
source("utils/plots/drug_response_pdx.R")

###########################################################
# Load in data
###########################################################

# pdx drug response data
xeva <- get_xeva("full")

# ------ get NK

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "komal"), ]

p_meta <- meta[meta$type == "PDX", ]

# load in arche scores
nk_pdxs_50k <- get_arche_scores("pdxs", "k50", p_meta) |> znorm()

# ------ get NT

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "tina"), ]

p_meta <- meta[meta$type == "PDX", ]

# load in arche scores
nt_pdxs_20k <- get_arche_scores("pdxs", "k20", p_meta) |> znorm()
nt_pdxs_50k <- get_arche_scores("pdxs", "k50", p_meta) |> znorm()

###########################################################
# Assign ARCHE scores
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
xeva_nk_pdxs_50k <- format_ARCHE(xeva, nk_pdxs_50k)
xeva_nt_pdxs_20k <- format_ARCHE(xeva, nt_pdxs_20k)
xeva_nt_pdxs_50k <- format_ARCHE(xeva, nt_pdxs_50k)

###########################################################
# Plot BAR
###########################################################

get_plot <- function(df, arche, drug) {

        subset_df <- df[df$drug == drug,]
        subset_df <- subset_df[
            complete.cases(subset_df[, c("BR_median", "BAR_median")]),
        ]
        p <- assess_ARCHE_TR(subset_df, arche, drug, "BAR_median", plot.indiv = TRUE)

}

get_plot(xeva_nk_pdxs_50k, "ARCHE4", "DATOPOTAMAB-CONTROL")
get_plot(xeva_nt_pdxs_50k, "ARCHE4", "SACITUZAMAB-GOVITECAN")
get_plot(xeva_nt_pdxs_50k, "ARCHE4", "AZD-8205")
get_plot(xeva_nt_pdxs_20k, "ARCHE4", "PACLITAXEL-15DAILY")
get_plot(xeva_nt_pdxs_20k, "ARCHE5", "CFI-400945-41.6WEEKLY")
get_plot(xeva_nt_pdxs_20k, "ARCHE5", "PACLITAXEL-15DAILY")