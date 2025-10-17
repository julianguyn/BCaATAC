# TODO: add in each step from NMF post-processing

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(data.table)
})

#import numpy as np
#data = np.load('6Rank_NNDSVD_Mixture.npy')
#np.savetxt('rank6Mixture.csv', data, delimiter=',')

###########################################################
# Load in data
###########################################################

# load in mixture from NMF
mixture <- as.data.frame(fread("/home/bioinf/bhklab/julia/projects/ATACseq/Signatures/results/ATAC_NMF_output/rank6Mixture.csv"))
rownames(mixture) <- paste0("Signature", 1:6)

# load in peak annotations
anno <- fread("peak_annotation.txt", header = F)$V1
colnames(mixture) <- anno

###########################################################
# Create BED file for background regions
###########################################################

split_strings <- strsplit(anno, ":")
background_BED <- as.data.frame(t(sapply(split_strings, function(x) unlist(x))))
colnames(background_BED) <- c("chrom", "chromStart", "chromEnd")
background_BED$chrom <- gsub("chr", "", background_BED$chrom)

write.table(background_BED, file = "beds/background.bed", quote = F, sep = "\t", col.names = T, row.names = F)

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
    filename <- paste0("beds/", signature, ".bed")
    write.table(BED, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
}