# temporary script, formating promoter file for HOMER to search models with UGT1A6 promoter accessibility

# load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

###########################################################
# Prepare promoter region data
###########################################################

# load in UGT1A6 promoter coordinates
coords <- read.csv("/cluster/home/julian/BCaATAC/data/UGT1A6.csv")

# turn into GRanges object
promoters <- GRanges(
    seqnames = "chr2",   
    ranges = IRanges(start = coords$Promoter.Start.Site,
                     end   = coords$Promoter.End.Site),
    gene = coords$Gene,
    transcript = coords$Ensembl
)

###########################################################
# Get BCa sample narrowpeak files
###########################################################

pdx_path <- "/cluster/projects/bhklab/projects/BCaATAC/BCa_ARCHE_Scoring/data/rawdata/PDXs/filtered_narrowpeak"
other_path <- "/cluster/projects/bhklab/projects/BCaATAC/BCa_ARCHE_Scoring/data/rawdata/CCLs/narrowpeaks"

pdx_files <- list.files(pdx_path, full.names = TRUE)
other_files <- list.files(other_path, full.names = TRUE)

peak_files <- c(pdx_files, other_files)
peak_files <- peak_files[-grep("sample_metadata", peak_files)]

peak_gr_list <- lapply(peak_files, function(f) {
    df <- fread(f, data.table = FALSE)
    colnames(df) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
    gr <- makeGRangesFromDataFrame(df, na.rm=TRUE)
    gr$sample <- basename(f)
    gr
})
names(peak_gr_list) <- basename(peak_files)


###########################################################
# Find overlaps
###########################################################

for (sample in names(peak_gr_list)) {
    sample_name <- sub("_peaks.*", "", sub("\\.", "", sample))
    print(paste("Working on", sample, "..."))
    peaks <- peak_gr_list[[sample]]
    hits <- findOverlaps(promoters, peaks)

    if (length(hits) == 0) {
        print("Found no hits")
    } else {
        overlap <- data.frame(
            sample = sample,
            gene = promoters$gene[queryHits(hits)],
            transcript = promoters$transcript[queryHits(hits)],
            promoter_start = start(promoters)[queryHits(hits)],
            promoter_end = end(promoters)[queryHits(hits)],
            peak_start = start(peaks)[subjectHits(hits)],
            peak_end = end(peaks)[subjectHits(hits)]
        )
        filename <- paste0("/cluster/projects/bhklab/projects/BCaATAC/Misc/data/results/", sample_name, ".rds")
        saveRDS(overlap, file = filename)
        print(paste("Saved hits to:", filename))
    }
}

###########################################################
# Merge overlap results
###########################################################

hits <- list.files("/cluster/projects/bhklab/projects/BCaATAC/Misc/data/results", full.names = TRUE)

compiled <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(compiled) <- c("sample", "gene", "transcript", "promoter_start", "promoter_end", "peak_start", "peak_end")

for (file in hits) {
    df <- readRDS(file)
    compiled <- rbind(compiled, df)
}

compiled$overlap_start <- pmax(compiled$promoter_start, compiled$peak_start)
compiled$overlap_end   <- pmin(compiled$promoter_end, compiled$peak_end)

# number of overlapping base pairs
compiled$overlap_bp <- pmax(0, compiled$overlap_end - compiled$overlap_start)

compiled$promoter_width <- compiled$promoter_end - compiled$promoter_start

# format dataframe to send to mitchell
compiled <- compiled[order(compiled$overlap_bp, decreasing = TRUE),]
compiled <- unique(compiled)

write.csv(compiled, file = "BCa_UGT1A6_overlap.csv", quote = FALSE, row.names = FALSE)