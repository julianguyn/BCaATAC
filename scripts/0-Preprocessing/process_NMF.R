# run scripts/1-Signatures/NMFcode/nmf.py first

# load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

source("utils/process_nmf.R")

###########################################################
# Load in data
###########################################################

# load in mixture from NMF
mixture <- as.data.frame(fread("data/procdata/NMF/rank6Mixture.csv"))
rownames(mixture) <- paste0("ARCHE", 1:6)

# get peak annotation
mat <- fread("data/rawdata/tcga/BCa_binary.2.matrix") |> as.data.frame()
anno <- paste(mat$V1, mat$V2, mat$V3, sep = ":")
write.table(anno, file = "data/procdata/ARCHEs/peak_annotation.txt")
colnames(mixture) <- anno

###########################################################
# Rank peaks per ARCHE
###########################################################

# order ARCHE peaks by weight and subset
#n <- 20000
n <- 50000

a1 <- top_peaks(mixture, "ARCHE1", n)
a2 <- top_peaks(mixture, "ARCHE2", n)
a3 <- top_peaks(mixture, "ARCHE3", n)
a4 <- top_peaks(mixture, "ARCHE4", n)
a5 <- top_peaks(mixture, "ARCHE5", n)
a6 <- top_peaks(mixture, "ARCHE6", n)

# plot changes in peak weights
plot_peakWeights(a1, "ARCHE1", n)
plot_peakWeights(a2, "ARCHE2", n)
plot_peakWeights(a3, "ARCHE3", n)
plot_peakWeights(a4, "ARCHE4", n)
plot_peakWeights(a5, "ARCHE5", n)
plot_peakWeights(a6, "ARCHE6", n)

###########################################################
# Create BED files for each ARCHE
###########################################################

# background BED file
createBED(anno, "Background")

# ARCHE bed files for top 20k sites
createBED(rownames(a1), "ARCHE1", n)
createBED(rownames(a2), "ARCHE2", n)
createBED(rownames(a3), "ARCHE3", n)
createBED(rownames(a4), "ARCHE4", n)
createBED(rownames(a5), "ARCHE5", n)
createBED(rownames(a6), "ARCHE6", n)

# ARCHE bed files for all
createBED(mixture, "ARCHE1", all = TRUE)
createBED(mixture, "ARCHE2", all = TRUE)
createBED(mixture, "ARCHE3", all = TRUE)
createBED(mixture, "ARCHE4", all = TRUE)
createBED(mixture, "ARCHE5", all = TRUE)
createBED(mixture, "ARCHE6", all = TRUE)

###########################################################
# Create HOMER peak files for each ARCHE
###########################################################

# background BED file for HOMER
createBEDforHOMER("Background")

# ARCHE HOMER peak files for top 20k sites
createBEDforHOMER("ARCHE1_20k")
createBEDforHOMER("ARCHE2_20k")
createBEDforHOMER("ARCHE3_20k")
createBEDforHOMER("ARCHE4_20k")
createBEDforHOMER("ARCHE5_20k")
createBEDforHOMER("ARCHE6_20k")

###########################################################
# Create Griffin site files for each ARCHE
###########################################################

createBEDforGriffin("ARCHE1_20k")
createBEDforGriffin("ARCHE2_20k")
createBEDforGriffin("ARCHE3_20k")
createBEDforGriffin("ARCHE4_20k")
createBEDforGriffin("ARCHE5_20k")
createBEDforGriffin("ARCHE6_20k")

###########################################################
# TODO: Create BEDs for BETA
###########################################################

# function to format bed file
formatBed <- function(bed) {
  colnames(bed) <- c("CHROM", "START", "END")
  bed$CHROM <- paste0("chr", bed$CHROM)
  bed$NAME <- paste(bed$CHROM, bed$START, bed$END, sep = ":")
  #bed$SCORE <- (nrow(bed):1) / 100
  bed$SCORE <- "."
  return(bed)
}