# process the TE bed files

#setwd("/cluster/projects/bhklab/projects/BCaATAC/peak-set-scoring/data/procdata/HERVs_LTRs")

count_regions <- data.frame(matrix(nrow=0, ncol=2))

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
    count_regions <- rbind(count_regions, data.frame(Family = dir, Regions = nrow(bed)))
}

saveRDS(count_regions, file = "count_regions.rds")


# individual bed files

for (file in list.files(recursive = TRUE, pattern = "*sorted")) {
    te <- read.table(file)
    colnames(te) <- c("seqnames", "start", "end")
    te$seqnames <- sub("chr", "", te$seqnames)
    te <- te[te$seqnames %in% c(1:22, "X", "Y"), ]
    filename <- paste0("all_beds/",sub(".*/", "", sub("\\.sorted", "", file)))
    write.table(te, file = filename, quote = FALSE, row.names = FALSE)
}

# stemness signatures

load("TE_stemness_signature.RData")

tp_bed <- data.frame(matrix(nrow=0, ncol=3))
pt_bed <- data.frame(matrix(nrow=0, ncol=3))

count_tp <- count_pt <- 0

for (file in list.files(recursive = TRUE, pattern = "*sorted")) {
    
    te <- sub(".*/", "", sub("\\.bed\\.sorted", "", file))

    if (te %in% te_tp) {
        te <- read.table(file)
        colnames(te) <- c("seqnames", "start", "end")
        te$seqnames <- sub("chr", "", te$seqnames)
        te <- te[te$seqnames %in% c(1:22, "X", "Y"), ]
        tp_bed <- rbind(tp_bed, te)
        count_tp <- count_tp + 1
    }
    if (te %in% te_pt) {
        te <- read.table(file)
        colnames(te) <- c("seqnames", "start", "end")
        te$seqnames <- sub("chr", "", te$seqnames)
        te <- te[te$seqnames %in% c(1:22, "X", "Y"), ]
        pt_bed <- rbind(pt_bed, te)
        count_pt <- count_pt + 1
    }

    write.table(tp_bed, file = "signature_bed/tp.bed", quote = FALSE, row.names = FALSE)
    write.table(pt_bed, file = "signature_bed/pt.bed", quote = FALSE, row.names = FALSE)
}

print(count_tp)
print(count_pt)