### Script to look at features of each signature (regions, annotate peaks, correlation, etc)

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


############################
### Correlate Signatures ###
############################ # run on server

setwd("/home/bioinf/bhklab/julia/projects/ATACseq")

library(data.table)
library(reshape2)
library(ggplot2)

# helper function to keep only upper triangle / lower triangle
# from: https://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# Get lower triangle of the correlation matrix
get_upper_tri<-function(cormat){
cormat[upper.tri(cormat)] <- NA
return(cormat)
}
# Get upper triangle of the correlation matrix
get_lower_tri <- function(cormat){
cormat[lower.tri(cormat)]<- NA
return(cormat)
}

# read in ATAC mixture
mixture <- fread("Signatures/results/ATAC_NMF_output/rank6Mixture.csv")
mixture <- t(mixture)
colnames(mixture) <- paste0("Signature", 1:6)

# compute correlation matrix
corr <- cor(mixture)

tri <- get_lower_tri(corr)
corr <- melt(tri)

# plot correlation matrix
png("Signatures/results/figures/ATAC_sig_corr.png", width = 6, height = 5, res = 600, units = "in", bg = "transparent")
ggplot(data = corr, aes(Var2, Var1, fill = value))+
    geom_tile(color = "transparent") + 
    scale_fill_gradient2(low = "#273469", high = "#6C2338", mid = "#F8E5EE", na.value = "transparent",
            midpoint = 0, limit = c(-1,1), space = "Lab", name="Correlation") +
    theme_void() + labs(y = "RNA-Seq Signatures\n", x = "\nATAC-Seq Signatures") +
    theme(axis.text.y = element_text(vjust = 0.5, hjust = 1),
          axis.title.y = element_text(angle = 90, margin = margin(t = 10)),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_text(margin = margin(t = 10)),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),) + coord_fixed()
dev.off()


# read in RNA mixture
mixture <- fread("Signatures/results/RNA_NMF_output/rank6Mixture.csv")
mixture <- t(mixture)
colnames(mixture) <- paste0("Signature", 1:6)

# compute correlation matrix
corr <- cor(mixture)

tri <- get_upper_tri(corr)
corr <- melt(tri)

# plot correlation matrix
png("Signatures/results/figures/RNA_sig_corr.png", width = 6, height = 5, res = 600, units = "in", bg = "transparent")
ggplot(data = corr, aes(Var2, Var1, fill = value))+
    geom_tile(color = "transparent") + 
    scale_fill_gradient2(low = "#273469", high = "#6C2338", mid = "#F8E5EE", na.value = "transparent",
            midpoint = 0, limit = c(-1,1), space = "Lab", name="Correlation") +
    theme_void() + labs(y = "RNA-Seq Signatures\n", x = "\nATAC-Seq Signatures") +
    theme(axis.text.y = element_text(vjust = 0.5, hjust = 1),
          axis.title.y = element_text(angle = 90, margin = margin(t = 10)),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_text(margin = margin(t = 10)),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),) + coord_fixed()
dev.off()