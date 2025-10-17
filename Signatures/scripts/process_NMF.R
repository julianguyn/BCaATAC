setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

#python stuff
#import numpy as np
#data = np.load('6Rank_NNDSVD_Mixture.npy')
#np.savetxt('rank6Mixture.csv', data, delimiter=',')
#data = np.load('6Rank_NNDSVD_Basis.npy')
#np.savetxt('rank6Basis.csv', data, delimiter=',')

###########################################################
# Load in data
###########################################################

# load in mixture from NMF
mixture <- as.data.frame(fread("Signatures/results/data/NMFoutputs/rank6Mixture.csv"))
rownames(mixture) <- paste0("ARCHE", 1:6)

# get peak annotation
mat <- fread("Signatures/data/BCa_binary.2.matrix") |> as.data.frame()
anno <- paste(mat$V1, mat$V2, mat$V3, sep = ":")
rm(mat)
write.table(anno, file = "Signatures/data/peak_annotation.txt")
colnames(mixture) <- anno

###########################################################
# Create BED file for background regions
###########################################################

split_strings <- strsplit(anno, ":")
background_BED <- as.data.frame(t(sapply(split_strings, function(x) unlist(x))))
colnames(background_BED) <- c("chrom", "chromStart", "chromEnd")
background_BED$chrom <- gsub("chr", "", background_BED$chrom)

write.table(background_BED, file = "Signatures/results/data/beds/background.bed", quote = F, sep = "\t", col.names = T, row.names = F)

###########################################################
# Rank peaks per ARCHE
###########################################################

# order ARCHE peaks by weight and subset
n <- 20000

a1 <- top_peaks(mixture, "ARCHE1", n)
a2 <- top_peaks(mixture, "ARCHE2", n)
a3 <- top_peaks(mixture, "ARCHE3", n)
a4 <- top_peaks(mixture, "ARCHE4", n)
a5 <- top_peaks(mixture, "ARCHE5", n)
a6 <- top_peaks(mixture, "ARCHE6", n)

# plot changes in peak weights
plot_peakWeights("ARCHE1", n)
plot_peakWeights("ARCHE2", n)
plot_peakWeights("ARCHE3", n)
plot_peakWeights("ARCHE4", n)
plot_peakWeights("ARCHE5", n)
plot_peakWeights("ARCHE6", n)

# ARCHE weight drops
arche1 <- c(947, 14241, 14351, 15991, 16493, 17459)

###########################################################
# Create BED file for each ARCHE
###########################################################

for (i in 1:nrow(mixture)) {
    tmp <- mixture[i,]

    # save signature name
    signature <- rownames(tmp)

    # keep only peaks with value
    tmp <- tmp[,-which(tmp==0)]

    # save peaks
    peaks <- colnames(tmp)

    # crate BED file
    split_strings <- strsplit(peaks, ":")
    BED <- as.data.frame(t(sapply(split_strings, function(x) unlist(x))))
    colnames(BED) <- c("chrom", "chromStart", "chromEnd")
    BED$chrom <- gsub("chr", "", BED$chrom)

    # save file
    filename <- paste0("Signatures/results/data/beds/", signature, ".bed")
    write.table(BED, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
}