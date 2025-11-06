# load libraries
suppressPackageStartupMessages({
    library(reshape2)
    library(data.table)
    library(ggplot2)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(GenomicRanges)
    library(dplyr)
    library(limma)
    library(DMRcate)
})

source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# get beta values
betas <- readRDS("data/procdata/TCGA/TCGA_betas.RDS")
rownames(betas) <- betas$CpG
betas$CpG <- NULL

# get 450k probe annotations
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# read in meta data file
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
meta$Sample.Name <- gsub("\\.", "-", meta$Sample.Name)
meta <- meta[match(colnames(betas), meta$Sample.Name), ]

###########################################################
# Remove ctrl and rs probes
###########################################################

betas <- betas[rownames(betas) %in% rownames(ann450k), ]

###########################################################
# Map CpGs to genomic coordinates
###########################################################

# helper function to find cpgs overlapping in arches
anno_cpgs <- function(arche) {

    # gr for arche
    bed <- fread(paste0("data/procdata/ARCHEs/beds/", arche, "_20k.bed"))
    gr <- GRanges(
        seqnames = paste0("chr", bed$chrom), 
        ranges = IRanges(bed$chromStart, bed$chromEnd)
    )

    # get overlap
    hits <- findOverlaps(anno_gr, gr)
    arche_cpgs <- ann450k[queryHits(hits), ] |> rownames()
    return(arche_cpgs)
}

# create gr of anno
anno_gr <- GRanges(
    seqnames = ann450k$chr,
    ranges = IRanges(start = ann450k$pos, end = ann450k$pos),
    probeID = rownames(ann450k)
)

# get cpgs in arche regions
arche1_cpgs <- anno_cpgs("ARCHE1")
arche2_cpgs <- anno_cpgs("ARCHE2")
arche3_cpgs <- anno_cpgs("ARCHE3")
arche4_cpgs <- anno_cpgs("ARCHE4")
arche5_cpgs <- anno_cpgs("ARCHE5")
arche6_cpgs <- anno_cpgs("ARCHE6")

###########################################################
# Distribution of beta values per ARCHE
###########################################################

# compile df of cpgs in arches
arche_cpgs <- data.frame(
    CpGs = c(arche1_cpgs, arche2_cpgs, arche3_cpgs, arche4_cpgs, arche5_cpgs, arche6_cpgs),
    ARCHE = c(
        rep("ARCHE1", length(arche1_cpgs)),
        rep("ARCHE2", length(arche2_cpgs)),
        rep("ARCHE3", length(arche3_cpgs)),
        rep("ARCHE4", length(arche4_cpgs)),
        rep("ARCHE5", length(arche5_cpgs)),
        rep("ARCHE6", length(arche6_cpgs))
    )
)

# create toPlot
toPlot <- data.frame(matrix(nrow=0, ncol=3))
for (arche in paste0("ARCHE", 1:6)) {
    samples <- meta$Sample.Name[meta$ARCHE == arche]

    # arche-specific cpg sites
    cpgs <- arche_cpgs$CpGs[arche_cpgs$ARCHE == arche]
    arche_bvals <- betas[rownames(betas) %in% cpgs, colnames(betas) %in% samples]

    # all cpgs for arche tumour
    all_bvals <- betas[, colnames(betas) %in% samples]
    
    bvals <- reshape2::melt(rbind(arche_bvals, all_bvals))
    bvals$ARCHE <- arche
    bvals$label <- c(
        rep("ARCHE\nspecific", ncol(arche_bvals)*nrow(arche_bvals)), 
        rep("All CpG\nsites", ncol(all_bvals)*nrow(all_bvals)))
    toPlot <- rbind(toPlot, bvals)
}

# categorize methylation status
toPlot <- toPlot %>%
    mutate(
        status = case_when(
            value <  0.2 ~ "Hypomethylated",
            value > 0.8 ~ "Hypermethylated",
            between(value, 0.2, 0.8) ~ "Intermediate",
            TRUE ~ NA_character_
        )
    )

# get proportions by methylation status
toPlot <- toPlot %>%
    filter(!is.na(status)) %>%
    group_by(ARCHE, label, status) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(ARCHE, label) %>%
    mutate(prop = n / sum(n))
toPlot$status <- factor(toPlot$status, levels = c("Hypomethylated", "Intermediate", "Hypermethylated"))


filename <- "data/results/figures/2-MolecularSigAnalysis/methylation/bval_arches.png"
png(filename, width = 6, height = 4, res = 600, units = "in")
ggplot(toPlot, aes(x = label, y = prop, fill = status)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ ARCHE) +
    geom_text(aes(label = paste0(round(prop * 100), "%")), position = position_fill(vjust = 0.5)) +
    scale_fill_manual(values = c(random_blue, "gray", random_lightblue)) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
    labs(x = "", y = "Proportion of CpG sites", fill = "Methylation Status")
dev.off()

###########################################################
# Convert beta values to M values
###########################################################

# add 1e-6 to avoid log 0
mvals <- log2((betas + 1e-6) / (1 - betas + 1e-6))

###########################################################
# Differentially methylated sites
###########################################################


###########################################################
# Differentially methylated regions
###########################################################