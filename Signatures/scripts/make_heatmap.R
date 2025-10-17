setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(reshape2)
    library(dplyr)
    library(RColorBrewer)
    library(ggh4x)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(ComplexHeatmap)
    library(GenomicRanges)
    library(wesanderson)
    library(rGREAT)
    library(RColorBrewer)
})

source("source/Signatures/plots.R")
source("source/palettes.R")

set.seed(123)

###########################################################
# Load in data
###########################################################

# read in meta data
meta <- read.csv("Signatures/data/TCGA_subtype_label.csv")

# load in matrix file from NMF
mat <- read.table("Signatures/results/data/ATAC_heatmap_rank6.png.order.matrix", header = T)

###########################################################
# Format signature data for plotting
###########################################################

# format for plotting and get subtype
mat$Signature <- paste0("ARCHE", 1:6)
mat$Signature <- factor(mat$Signature, levels = paste0("ARCHE",6:1))
mat <- reshape2::melt(mat)
mat$subtype <- meta$bca_subtype[match(gsub("X", "", mat$variable), gsub("-", "\\.", meta$File.Name))]

# save signature assignment
mat$signature_assign <- ""
for (sample in mat$variable) {
    tmp <- mat[mat$variable == sample,]
    mat[mat$variable == sample,]$signature_assign <- as.character(tmp[which.max(tmp$value),]$Signature)
}

###########################################################
# Plot ATAC-Signature heatmap
###########################################################

plot_ARCHE_heatmap(mat)

###########################################################
# Load in ARCHE BEDs
###########################################################

sig1 <- fread("Signatures/results/data/beds/Signature1.bed")
sig2 <- fread("Signatures/results/data/beds/Signature2.bed")
sig3 <- fread("Signatures/results/data/beds/Signature3.bed")
sig4 <- fread("Signatures/results/data/beds/Signature4.bed")
sig5 <- fread("Signatures/results/data/beds/Signature5.bed")
sig6 <- fread("Signatures/results/data/beds/Signature6.bed")

###########################################################
# Compute number and size of peak regions
###########################################################

# compute total open peak region
sig1$diff <- sig1$chromEnd - sig1$chromStart
sig2$diff <- sig2$chromEnd - sig2$chromStart
sig3$diff <- sig3$chromEnd - sig3$chromStart
sig4$diff <- sig4$chromEnd - sig4$chromStart
sig5$diff <- sig5$chromEnd - sig5$chromStart
sig6$diff <- sig6$chromEnd - sig6$chromStart

# get peak information
num_windows <- c(nrow(sig1), nrow(sig2), nrow(sig3), nrow(sig4), nrow(sig5), nrow(sig6))
sum_peaks <- c(sum(sig1$diff), sum(sig2$diff), sum(sig3$diff), sum(sig4$diff), sum(sig5$diff), sum(sig6$diff))

df <- data.frame(Signature = paste0("Signature", 1:6),
                num_windows = num_windows, sum_peaks = sum_peaks)

###########################################################
# Compute number of overlapping regions
###########################################################

# create GRanges
gr1 <- GRanges(seqnames = paste0("chr", sig1$chrom), ranges = IRanges(sig1$chromStart, sig1$chromEnd))
gr2 <- GRanges(seqnames = paste0("chr", sig2$chrom), ranges = IRanges(sig2$chromStart, sig2$chromEnd))
gr3 <- GRanges(seqnames = paste0("chr", sig3$chrom), ranges = IRanges(sig3$chromStart, sig3$chromEnd))
gr4 <- GRanges(seqnames = paste0("chr", sig4$chrom), ranges = IRanges(sig4$chromStart, sig4$chromEnd))
gr5 <- GRanges(seqnames = paste0("chr", sig5$chrom), ranges = IRanges(sig5$chromStart, sig5$chromEnd))
gr6 <- GRanges(seqnames = paste0("chr", sig6$chrom), ranges = IRanges(sig6$chromStart, sig6$chromEnd))

# create peak list for Upset plot
peak_list <- list(
    ARCHE1 = gr1, 
    ARCHE2 = gr2, 
    ARCHE3 = gr3, 
    ARCHE4 = gr4, 
    ARCHE5 = gr5, 
    ARCHE6 = gr6
)

# plot UPSET plot of overlapping peaks
m = make_comb_mat(peak_list)
m = m[comb_size(m) > 2000000]
plot_ATAC_Upset(m)

###########################################################
# Annotate ARCHE peak sets
###########################################################

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# annotate peaks
anno1 <- annotatePeak(gr1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
anno2 <- annotatePeak(gr2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
anno3 <- annotatePeak(gr3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
anno4 <- annotatePeak(gr4, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
anno5 <- annotatePeak(gr5, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
anno6 <- annotatePeak(gr6, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat

# add ARCHE labels
anno1$ARCHE <- "ARCHE1"
anno2$ARCHE <- "ARCHE2"
anno3$ARCHE <- "ARCHE3"
anno4$ARCHE <- "ARCHE4"
anno5$ARCHE <- "ARCHE5"
anno6$ARCHE <- "ARCHE6"

toPlot <- rbind(anno1, anno2, anno3, anno4, anno5, anno6)

png("Signatures/results/figures/peakAnno.png", width = 7, height = 5, res = 600, units = "in")
ggplot(toPlot, aes(fill=Feature, y=Frequency, x=Signature)) + 
    geom_bar(position="fill", stat="identity", color = "black", size = 0.3) +
    scale_fill_manual(values = brewer.pal(11, name = "Paired")) +
    theme_minimal() + labs(y = "Percentage (%)")
dev.off()