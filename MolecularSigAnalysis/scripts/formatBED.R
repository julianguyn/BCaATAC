# load libraries
suppressPackageStartupMessages({
  library(data.table)
})

# script to format BED files for BETA


###########################################################
# Load in data
###########################################################

# load in counts matrix
sig1 <- fread("Signatures/Signature1.bed")
sig2 <- fread("Signatures/Signature2.bed")
sig3 <- fread("Signatures/Signature3.bed")
sig4 <- fread("Signatures/Signature4.bed")
sig5 <- fread("Signatures/Signature5.bed")
sig6 <- fread("Signatures/Signature6.bed")


###########################################################
# Format peak file for BETA
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

bed1 <- formatBed(sig1)
bed2 <- formatBed(sig2)
bed3 <- formatBed(sig3)
bed4 <- formatBed(sig4)
bed5 <- formatBed(sig5)
bed6 <- formatBed(sig6)


###########################################################
# Write peak bed files
###########################################################

write.table(bed1, "MolecularSigAnalysis/data/bed1.bed",row.names = F,col.names = F, sep="\t", quote=FALSE) 
write.table(bed2, "MolecularSigAnalysis/data/bed2.bed",row.names = F,col.names = F, sep="\t", quote=FALSE) 
write.table(bed3, "MolecularSigAnalysis/data/bed3.bed",row.names = F,col.names = F, sep="\t", quote=FALSE) 
write.table(bed4, "MolecularSigAnalysis/data/bed4.bed",row.names = F,col.names = F, sep="\t", quote=FALSE) 
write.table(bed5, "MolecularSigAnalysis/data/bed5.bed",row.names = F,col.names = F, sep="\t", quote=FALSE) 
write.table(bed6, "MolecularSigAnalysis/data/bed6.bed",row.names = F,col.names = F, sep="\t", quote=FALSE) 