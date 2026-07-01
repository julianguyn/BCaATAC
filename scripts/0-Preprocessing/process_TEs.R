# process the TE bed files

setwd("/cluster/projects/bhklab/projects/BCaATAC/peak-set-scoring/data/procdata/HERVs_LTRs")

for (dir in list.files()) {

    bed <- data.frame(matrix(nrow=0, ncol=3))

    for (file in list.files(dir)) {
        te <- read.table(paste0(dir, "/", file))
        colnames(te) <- c("seqnames", "start", "end")
        te$seqnames <- sub("chr", "", te$seqnames)
        te <- te[te$seqnames %in% c(1:22, "X", "Y"), ]
        bed <- rbind(bed, te)
    }
    write.table(bed, file = paste0("beds/", dir, ".bed"), quote = FALSE, row.names = FALSE)
}