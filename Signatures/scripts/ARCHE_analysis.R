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
    library(rGREAT)
    library(RColorBrewer)
})

source("source/Signatures/plots.R")
source("source/Signatures/helper.R")
source("source/palettes.R")

set.seed(123)

# ******************************
# specify which BED files to use
# ******************************
analysis = c("20k", "all") 

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
# Load in BEDs
###########################################################

bg <- fread(paste0("Signatures/results/data/beds/Background.bed"))

sig1 <- fread(paste0("Signatures/results/data/beds/ARCHE1_", analysis, ".bed"))
sig2 <- fread(paste0("Signatures/results/data/beds/ARCHE2_", analysis, ".bed"))
sig3 <- fread(paste0("Signatures/results/data/beds/ARCHE3_", analysis, ".bed"))
sig4 <- fread(paste0("Signatures/results/data/beds/ARCHE4_", analysis, ".bed"))
sig5 <- fread(paste0("Signatures/results/data/beds/ARCHE5_", analysis, ".bed"))
sig6 <- fread(paste0("Signatures/results/data/beds/ARCHE6_", analysis, ".bed"))

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

# plot peak info
df <- data.frame(ARCHE = paste0("ARCHE", 1:6),
                num_windows = num_windows, sum_peaks = sum_peaks)
plot_ARCHE_peakInfo(df)

###########################################################
# Compute number of overlapping regions
###########################################################

# create GRanges
grb <- GRanges(seqnames = paste0("chr", bg$chrom), ranges = IRanges(bg$chromStart, bg$chromEnd))

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
plot_ATAC_Upset(m, analysis)

###########################################################
# Annotate ARCHE peak sets
###########################################################

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# annotate peaks
anno1 <- annotateARCHE(gr1, "ARCHE1")
anno2 <- annotateARCHE(gr2, "ARCHE2")
anno3 <- annotateARCHE(gr3, "ARCHE3")
anno4 <- annotateARCHE(gr4, "ARCHE4")
anno5 <- annotateARCHE(gr5, "ARCHE5")
anno6 <- annotateARCHE(gr6, "ARCHE6")
annob <- annotateARCHE(grb, "Background")

# plot peakAnno results
toPlot <- rbind(anno1, anno2, anno3, anno4, anno5, anno6, annob)
plot_annotatePeak(toPlot, analysis)

###########################################################
# GREAT analysis
###########################################################

# run GREAT
great1 <- runGREAT(gr1, "ARCHE1", analysis)
great2 <- runGREAT(gr2, "ARCHE2", analysis)
great3 <- runGREAT(gr3, "ARCHE3", analysis)
great4 <- runGREAT(gr4, "ARCHE4", analysis)
great5 <- runGREAT(gr5, "ARCHE5", analysis)
great6 <- runGREAT(gr6, "ARCHE6", analysis)


# plot results
plot_GREAT <- function(great) {

    # keep only top 30
    toPlot <- df[1:30,]
    toPlot <- toPlot[order(toPlot$Adjp_BH, decreasing = TRUE),]
    toPlot$GO.Name <- factor(toPlot$GO.Name, levels=unique(toPlot$GO.Name))

    p <- ggplot(toPlot, aes(x = GO.Name, y = -log(Adjp_BH), size = Total_Genes_Annotated, color = Label), shape = 19) +
        geom_point() + #scale_size(range = c(0.1, 3)) +
        coord_cartesian(clip = "off") + coord_flip()  +
        guides(shape = guide_legend(ncol = 1), color = guide_legend(override.aes=list(shape=19, size = 4))) +
        #guides(shape = guide_legend(override.aes = list(size = 1)), color = guide_legend(override.aes = list(shape = 19, size = 1))) +
        scale_color_manual(values = c("#574B60", "#CBC9AD", "#5B96AF"), labels = c("BP", "CC", "MF")) +
        theme_classic() + 
        theme(legend.key.size = unit(0.5, 'cm'), text = element_text(size = 10),
            legend.position = c(0.76, 0.2), #plot.margin=unit(c(0.6,0,0,1),"cm"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + 
        labs(x = "GO Term", y = "-log(FDR)", size = "N", color = "Ontology")
    return(p)
}