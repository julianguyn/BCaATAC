# script to format tracks from Chip-Atlas for chromVar

library(data.table)

# load in and format file
bed <- fread("data/rawdata/tracks/Oth.Brs.05.NFE2L2.AllCell.bed", skip = "track", header = FALSE, sep = "\t")
bed <- bed[,c(1:3)]
colnames(bed) <- c("seqnames", "start", "end")
bed$seqnames <- sub("chr", "", bed$seqnames)

# write file
filename <- paste0("data/rawdata/tracks/NRF2_TF_sites.bed")
write.table(bed, file = filename, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
