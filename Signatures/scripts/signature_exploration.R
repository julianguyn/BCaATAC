### Script to look at features of each signature (regions, annotate peaks, etc)

setwd("C:/Users/julia/Documents/BCaATAC")

suppressMessages(library(rGREAT))
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ComplexHeatmap)
library(GenomicRanges)
library(data.table)
library(wesanderson)
library(ggplot2)

set.seed(123)

# read in signature bed files
sig1 <- fread("Signatures/results/data/beds/Signature1.bed")
sig2 <- fread("Signatures/results/data/beds/Signature2.bed")
sig3 <- fread("Signatures/results/data/beds/Signature3.bed")
sig4 <- fread("Signatures/results/data/beds/Signature4.bed")
sig5 <- fread("Signatures/results/data/beds/Signature5.bed")
sig6 <- fread("Signatures/results/data/beds/Signature6.bed")


###############################################
### Compute number and size of peak regions ###
###############################################

# compute total open peak region
sig1$diff <- sig1$chromEnd - sig1$chromStart
sig2$diff <- sig2$chromEnd - sig2$chromStart
sig3$diff <- sig3$chromEnd - sig3$chromStart
sig4$diff <- sig4$chromEnd - sig4$chromStart
sig5$diff <- sig5$chromEnd - sig5$chromStart
sig6$diff <- sig6$chromEnd - sig6$chromStart

# get peak information
num_peaks <- c(nrow(sig1), nrow(sig2), nrow(sig3), nrow(sig4), nrow(sig5), nrow(sig6))
sum_peaks <- c(sum(sig1$diff), sum(sig2$diff), sum(sig3$diff), sum(sig4$diff), sum(sig5$diff), sum(sig6$diff))
avg_peak <- c(mean(sig1$diff), mean(sig2$diff), mean(sig3$diff), mean(sig4$diff), mean(sig5$diff), mean(sig6$diff))

df <- data.frame(Signature = paste0("Signature", 1:6),
                num_peaks = num_peaks, sum_peaks = sum_peaks, avg_peak = avg_peak)


#############################################
### Compute number of overlapping regions ###
#############################################

# convert files to GRanges object
gr1 <- GRanges(seqnames = sig1$chrom, ranges = IRanges(sig1$chromStart, sig1$chromEnd))
gr2 <- GRanges(seqnames = sig2$chrom, ranges = IRanges(sig2$chromStart, sig2$chromEnd))
gr3 <- GRanges(seqnames = sig3$chrom, ranges = IRanges(sig3$chromStart, sig3$chromEnd))
gr4 <- GRanges(seqnames = sig4$chrom, ranges = IRanges(sig4$chromStart, sig4$chromEnd))
gr5 <- GRanges(seqnames = sig5$chrom, ranges = IRanges(sig5$chromStart, sig5$chromEnd))
gr6 <- GRanges(seqnames = sig6$chrom, ranges = IRanges(sig6$chromStart, sig6$chromEnd))

# create peak list for Upset plot
peak_list <- list(Signature1 = gr1, Signature2 = gr2, Signature3 = gr3, 
                  Signature4 = gr4, Signature5 = gr5, Signature6 = gr6)


m = make_comb_mat(peak_list)
m = m[comb_size(m) > 2000000]

png("Signatures/results/ATAC_Upset.png", width = 11, height = 5, res = 600, units = "in")
UpSet(m, set_order = c(paste0("Signature", 1:6)), comb_order = order(-comb_size(m)),
      bg_col = "#F5F0EC", comb_col = "#463239", 
      right_annotation = upset_right_annotation(m, gp = gpar(fill = "#8F8073")))
dev.off()


####################################
### Annotate Signature Peak Sets ###
####################################

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

gr1 <- GRanges(seqnames = paste0("chr", sig1$chrom), ranges = IRanges(sig1$chromStart, sig1$chromEnd))
gr2 <- GRanges(seqnames = paste0("chr", sig2$chrom), ranges = IRanges(sig2$chromStart, sig2$chromEnd))
gr3 <- GRanges(seqnames = paste0("chr", sig3$chrom), ranges = IRanges(sig3$chromStart, sig3$chromEnd))
gr4 <- GRanges(seqnames = paste0("chr", sig4$chrom), ranges = IRanges(sig4$chromStart, sig4$chromEnd))
gr5 <- GRanges(seqnames = paste0("chr", sig5$chrom), ranges = IRanges(sig5$chromStart, sig5$chromEnd))
gr6 <- GRanges(seqnames = paste0("chr", sig6$chrom), ranges = IRanges(sig6$chromStart, sig6$chromEnd))

anno1 <- annotatePeak(gr1, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
anno2 <- annotatePeak(gr2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
anno3 <- annotatePeak(gr3, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
anno4 <- annotatePeak(gr4, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
anno5 <- annotatePeak(gr5, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
anno6 <- annotatePeak(gr6, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat

anno1$Signature <- "Signature1"
anno2$Signature <- "Signature2"
anno3$Signature <- "Signature3"
anno4$Signature <- "Signature4"
anno5$Signature <- "Signature5"
anno6$Signature <- "Signature6"

toPlot <- rbind(anno1, anno2, anno3, anno4, anno5, anno6)

library(RColorBrewer)

png("Signatures/results/figures/peakAnno.png", width = 7, height = 5, res = 600, units = "in")
ggplot(toPlot, aes(fill=Feature, y=Frequency, x=Signature)) + 
    geom_bar(position="fill", stat="identity", color = "black", size = 0.3) +
    scale_fill_manual(values = brewer.pal(11, name = "Paired")) +
    theme_minimal() + labs(y = "Percentage (%)")
dev.off()