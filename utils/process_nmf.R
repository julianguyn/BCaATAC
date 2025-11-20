#' Order ARCHE peaks by weight and subset for top n peaks
#' 
#' @param mixture dataframe. From NMF
#' @param arche string. ARCHE name
#' @param num_windows int. Number of windows to subset for
#' @return dataframe of top n peaks in ARCHE
#' 
top_peaks <- function(mixture, arche, num_windows) {

    peaks <- mixture[rownames(mixture) == arche,]
    peaks <- peaks[,-which(peaks==0)] |> t() |> as.data.frame()
    peaks <- peaks[order(peaks[[arche]], decreasing = T),,drop=F]

    # subset top windows
    peaks <- peaks[1:num_windows,,drop=F]
    peaks$rank <- 1:nrow(peaks)

    return(peaks)
}

#' Find locations of plateaus
#' 
#' @param df data.frame. From top_peaks()
#' @return dataframe of start positions of plateaus of >100 windows
#' 
find_delta0 <- function(df) {

  # get difference between sequential weights
  colnames(df)[1] <- "Weight"
  df <- df %>%
    mutate(diff = lead(Weight) - Weight)

  # remove plateaus
  df <- df[df$diff < 0,]

  #get difference between ranks
  df <- df %>%
    mutate(diff_rank = lead(rank) - rank)
  
  # get start positions of plateaus of >100 windows
  df <- df[which(df$diff_rank > 100),]

  return(df)
}

#' ARCHE peak weight thresholds
#' 
#' Plot changes in peak weights from NMF for top n windows
#' @param peaks dataframe. Output of top_peaks()
#' @param arche string. ARCHE name
#' @param num_windows int. Number of windows kept
#' 
plot_peakWeights <- function(peaks, arche, num_windows) {

    n <- paste0(as.character(num_windows/1000), "k")

    # get start positions of plateaus
    diff <- find_delta0(peaks)

    filename <- paste0("data/results/figures/0-Preprocessing/",arche,"_", n, ".png")
    png(filename, width = 5, height = 5, res = 600, units = "in")
    print(
        ggplot(peaks, aes(x = rank, y = .data[[arche]])) +
        geom_line() +
        geom_point() +
        geom_vline(xintercept = diff$rank, linetype = "dashed") +
        geom_vline(xintercept = 10000) +
        theme_classic() +
        labs(x = "Peak Window Ranked by Weight", y = paste("Peak Window Weight for", arche))
    )
    dev.off()
}

#' Create BED file 
#' 
#' @param peaks string. Vector of coordinates chr:start:end
#' @param filename string. Filename for BED
#' 
createBED <- function(peaks, filename, n = 0, all = FALSE) {

  n <- ifelse(
    all == TRUE,
    "all",
    paste0(as.character(n/1000), "k")
  )

  # if all == TRUE
  if (all == TRUE) {
    arche <- gsub("_all", "", filename)
    peaks <- peaks[rownames(peaks) == arche,]
    peaks <- peaks[,-which(peaks==0)]
    peaks <- colnames(peaks)
  }

  # create BED
  split_strings <- strsplit(peaks, ":")
  BED <- as.data.frame(t(sapply(split_strings, function(x) unlist(x))))
  colnames(BED) <- c("chrom", "chromStart", "chromEnd")
  BED$chrom <- gsub("chr", "", BED$chrom)

  # save file
  filename <- paste0("data/procdata/ARCHEs/beds/", filename, "_", n, ".bed")
  write.table(BED, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
}

#' Create HOMER peak file 
#' 
#' @param filename string. Same was used in createBED()
#' 
createBEDforHOMER <- function(filename) {

  bed <- fread(paste0("data/procdata/ARCHEs/beds/", filename, ".bed"), data.table=F)
  peak <- data.frame(
    Chr = paste0("chr", bed$chrom),
    Start = bed$chromStart,
    End = bed$chromEnd,
    PeakID = paste(bed$chrom, bed$chromStart, bed$chromEnd, sep = ":"),
    Skip = "pleasework",
    Strand = "."
  )
  # save file
  filename <- paste0("data/procdata/ARCHEs/HOMERbeds/", filename, ".bed")
  write.table(peak, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
}