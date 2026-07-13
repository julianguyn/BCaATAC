# check for stemness TE signatures

# analyze individual groups of TEs

suppressPackageStartupMessages({
    library(readxl)
    library(tidyverse)
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    library(reshape2)
    library(ggnewscale)
})

source("utils/get_data.R")
source("utils/palettes.R")

set.seed(123)

## --------------------------------------------------------
# Preprocessing for H4H

# load in stemness signatures
te_tp <- read_excel("data/rawdata/TEs/cd-23-0331_supplementary_table_s2_suppst2.xlsx", sheet = 1, col_names = TRUE, skip = 1)
te_pt <- read_excel("data/rawdata/TEs/cd-23-0331_supplementary_table_s3_suppst3.xlsx", sheet = 1, col_names = TRUE, skip = 1)

te_tp <- te_tp$'TE family'[te_tp$Enrichment == "Mammary"] #56
te_pt <- te_pt$'TE family'[te_pt$Enrichment == "PSCs_Mammary"] #194

save(te_tp, te_pt, file = "data/procdata/TCGA/TE_stemness_signature.RData")

# from here, run 0-Preprocessing/process_TEs.R

###########################################################
# Load in data
###########################################################

# load in metadata
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
meta$Sample.Name <- gsub("\\.", "-", meta$Sample.Name)

# get signature scores from NMF
mat <- get_arche_tcga()
mat_order <- mat$variable
mat$variable <- meta$Sample.Name[match(mat$variable, meta$ATAC.Seq.File.Name)]

# load in zscores
tes <- read.table("data/rawdata/TCGA_TEs/tcga_TE_stemness.Zscore.txt")

###########################################################
# Handle duplicates and format data
###########################################################

# rename duplicate
mat[mat$variable == "TCGA-A2-A0T4",]$variable <- rep(paste0("TCGA-A2-A0T4-", 1:2), each = 6)
mapping <- unique(data.frame(mat = mat_order, new = mat$variable))
colnames(tes) <- mapping$new[match(sub("X", "", colnames(tes)), mapping$mat)]

###########################################################
# Get variables
###########################################################

tcga_anno <- unique(mat[,c(2,4,5)])

###########################################################
# Create toPlot
###########################################################

tes$Signature <- rownames(tes)
toPlot <- reshape2::melt(tes)
toPlot$ARCHE <- mat$signature_assign[match(toPlot$variable, mat$variable)]
toPlot$Subtype <- mat$subtype[match(toPlot$variable, mat$variable)]

###########################################################
# Plot
###########################################################

p <- ggplot(toPlot, aes(x = ARCHE, y = value)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(aes(fill = Signature), outlier.shape = NA,
                 position = position_dodge(width = 0.75)) +
    scale_fill_manual("TE Signature", 
        values = c("pt" = "#DFDFDF", "tp" = "#949396"),
        labels = c("PSC-enriched", "Tissue-enriched")
    ) +
    new_scale_fill() +
    geom_jitter(aes(fill = Subtype, group = Signature), 
                shape = 21, colour = "black", size = 2,
                position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75)) +
    scale_fill_manual("Subtype", values = subtype_pal) +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    labs(y = "Zscore") +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        legend.key.size = unit(0.5, 'cm')
    )

filename <- "data/results/figures/1-Signatures/TE_stemness.png"
ggsave(filename, p, width = 6, height = 3)
