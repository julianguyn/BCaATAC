# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(grid)
    library(gridExtra)
})


source("utils/palettes.R")
source("utils/get_data.R")
source("utils/mappings.R")
source("utils/plots/drug_response_pdx_indivplots.R")
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

p_meta_t <- meta_t[meta_t$type == "PDX", ]
p_meta_k <- meta_k[meta_k$type == "PDX", ]

###########################################################
# Load in data
###########################################################

# read in RNA
rna <- read.csv("data/rawdata/pdx/gene_tpm_normalized_matrix.csv")
rownames(rna) <- rna$X
colnames(rna) <- sub("^S", "", colnames(rna))
colnames(rna) <- map_pdx(colnames(rna))
rna <- rna[,colnames(rna) %in% p_meta_t$sampleid]

# load in PDX drug response (from "4-DrugResponse/ARCHE_PDXs.R")
xeva <- get_xeva("full")

# load in PDX subtyping (from 3-DataExploration/PDXs_rna.R)
load("data/results/data/3-DataExploration/pdxs_subtyping_scores.RData")
# pdx_pam50

###########################################################
# Add in subtype information
###########################################################

# helper function
format_ARCHE <- function(df, meta) {
    df$Subtype <- NA
    for (i in 1:nrow(df)) {
        sample <- df$patient.id[i]
        if (sample %in% meta$sampleid) {
            df$Subtype[i] <- meta$subtype[meta$sampleid == sample]
        } else {
            df$Subtype[i] <- "Not Available"
        }
    }
    # manual mapping
    df$Subtype[df$PDX_ID == "REF032"] <- "ER"
    df$Subtype[df$PDX_ID == "REF034"] <- "ER"
    df$Subtype[df$PDX_ID == "REF036"] <- "ER"
    df$Subtype[df$PDX_ID == "REF038"] <- "TNBC"
    df$Subtype[df$PDX_ID == "NOTCH01"] <- "TNBC"
    df$Subtype[df$Subtype == "unknown"] <- "Not Available"
    return(df)
}

xeva_t <- format_ARCHE(xeva, p_meta_t)
xeva_k <- format_ARCHE(xeva, p_meta_k)

###########################################################
# Plot paclitaxel TR associations
###########################################################

plot_paclitaxel_TR <- function(xeva, label) {
    subset_df <- xeva[xeva$drug == "PACLITAXEL-15DAILY",]
    subset_df$Basal <- pdx_pam50$Basal[match(subset_df$PDX_ID, rownames(pdx_pam50))]
    subset_df <- subset_df[
        complete.cases(subset_df[, c("BR_median", "BAR_median")]),
    ]
    assess_ARCHE_TR(subset_df, "Basal", "PACLITAXEL-15DAILY", "BR_median", label = label, plot.indiv = T)
    assess_ARCHE_TR(subset_df, "Basal", "PACLITAXEL-15DAILY", "BAR_median", label = label, plot.indiv = T)
}

plot_paclitaxel_TR(xeva_t, "xeva_t")
plot_paclitaxel_TR(xeva_k, "xeva_k")

###########################################################
# Plot ADC TOPIi TR associations
###########################################################
