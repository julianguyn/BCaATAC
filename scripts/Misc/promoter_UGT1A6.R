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

path <- "/cluster/projects/bhklab/rawdata/BCaATAC/new_files"
pdx_files <- list.files(path, full.names = TRUE, recursive = TRUE, pattern = ".*peaks.filtered.narrowPeak")
cell_files <- list.files(path, full.names = TRUE, recursive = TRUE, pattern = ".*peaks.filtered.merged.narrowPeak")

peak_files <- c(pdx_files, cell_files)

peak_gr_list <- lapply(peak_files, function(f) {
    df <- fread(f, data.table = FALSE)
    df <- df[,c(1:3)]
    colnames(df) <- c("chrom", "chromStart", "chromEnd")
    gr <- makeGRangesFromDataFrame(df, na.rm = TRUE)
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

hits <- list.files("/cluster/projects/bhklab/projects/BCaATAC/Misc/data/results", full.names = TRUE, pattern = ".*.rds")

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

write.csv(compiled, file = "/cluster/projects/bhklab/projects/BCaATAC/Misc/data/results/BCa_UGT1A6_overlap.csv", quote = FALSE, row.names = FALSE)


# send to Mitchell:
compiled[compiled$overlap_bp == 1200, c(1:3)]
                                            sample   gene        transcript
73    66684_P1_S6_peaks.filtered.merged.narrowPeak UGT1A6   ENSG00000167165
74    66684_P1_S6_peaks.filtered.merged.narrowPeak UGT1A6 ENST00000373424.5
77    66684_P1_S6_peaks.filtered.merged.narrowPeak UGT1A6 ENST00000441351.1
111       BT20_1_S1_L001_peaks.filtered.narrowPeak UGT1A6   ENSG00000167165
117       BT20_1_S1_L001_peaks.filtered.narrowPeak UGT1A6 ENST00000373424.5
122       BT20_1_S1_L001_peaks.filtered.narrowPeak UGT1A6 ENST00000441351.1
125     BT20_1_S1_peaks.filtered.merged.narrowPeak UGT1A6   ENSG00000167165
127     BT20_1_S1_peaks.filtered.merged.narrowPeak UGT1A6 ENST00000373424.5
129     BT20_1_S1_peaks.filtered.merged.narrowPeak UGT1A6 ENST00000441351.1
262   HCC1806_2_S21_L002_peaks.filtered.narrowPeak UGT1A6   ENSG00000167165
268   HCC1806_2_S21_L002_peaks.filtered.narrowPeak UGT1A6 ENST00000373424.5
274   HCC1806_2_S21_L002_peaks.filtered.narrowPeak UGT1A6 ENST00000441351.1
279 HCC1806_2_S21_peaks.filtered.merged.narrowPeak UGT1A6   ENSG00000167165
282 HCC1806_2_S21_peaks.filtered.merged.narrowPeak UGT1A6 ENST00000373424.5
285 HCC1806_2_S21_peaks.filtered.merged.narrowPeak UGT1A6 ENST00000441351.1
310    HCC1954_1_S9_L001_peaks.filtered.narrowPeak UGT1A6   ENSG00000167165
312    HCC1954_1_S9_L001_peaks.filtered.narrowPeak UGT1A6 ENST00000373424.5
317    HCC1954_1_S9_L001_peaks.filtered.narrowPeak UGT1A6 ENST00000441351.1
319  HCC1954_1_S9_peaks.filtered.merged.narrowPeak UGT1A6   ENSG00000167165
320  HCC1954_1_S9_peaks.filtered.merged.narrowPeak UGT1A6 ENST00000373424.5
322  HCC1954_1_S9_peaks.filtered.merged.narrowPeak UGT1A6 ENST00000441351.1
507    REF024_S17_peaks.filtered.merged.narrowPeak UGT1A6 ENST00000406651.1