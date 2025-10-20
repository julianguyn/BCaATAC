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

#' Create BED file 
#' 
#' @param peaks string. Vector of coordinates chr:start:end
#' @param filename string. Filename for BED
#' 
createBED <- function(peaks, filename, all = FALSE) {

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
  filename <- paste0("Signatures/results/data/beds/", filename, ".bed")
  write.table(BED, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
}

#' Create HOMER peak file 
#' 
#' @param filename string. Same was used in createBED()
#' 
createBEDforHOMER <- function(filename) {

  bed <- fread(paste0("Signatures/results/data/beds/", filename, ".bed"), data.table=F)
  peak <- data.frame(
    Chr = paste0("chr", bed$chrom),
    Start = bed$chromStart,
    End = bed$chromEnd,
    PeakID = paste(bed$chrom, bed$chromStart, bed$chromEnd, sep = ":"),
    Skip = "pleasework",
    Strand = "."
  )
  # save file
  filename <- paste0("Signatures/results/data/HOMERbeds/", filename, ".bed")
  write.table(peak, file = filename, quote = F, sep = "\t", col.names = T, row.names = F)
}

#' Annotate peaks in ARCHE
#' 
#' @param gr GRanges object. GRanges of ARCHE peaks
#' @param arche string. ARCHE label
#' @return anno annotatePeak result
#' 
annotateARCHE <- function(gr, arche) {

    anno <- annotatePeak(gr, tssRegion=c(-250, 250), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
    anno$ARCHE <- arche
    return(anno)
}

#' Run GREAT on ARCHE regions
#' 
#' @param gr GRanges object. GRanges of ARCHE peaks
#' @param arche string. ARCHE label
#' @param analysis string. 20k or all for label
#' 
runGREAT <- function(gr, arche, analysis) {
    
    message("****Job submitted")
    job = submitGreatJob(gr, bg, species = "hg38", genome = "hg38", help = FALSE)
    message("****Job completed")
    tbl = getEnrichmentTables(job)

    # save results
    mf <- as.data.frame(tbl[1])
    bp <- as.data.frame(tbl[2])
    cc <- as.data.frame(tbl[3])

    # save column names
    cols <- gsub("GO.Molecular.Function.", "", colnames(mf))
    colnames(mf) <- colnames(bp) <- colnames(cc) <- cols

    mf$Label <- "Molecular Feature"
    bp$Label <- "Biological Process"
    cc$Label <- "Cellular Component"

    # combine results and filter
    res <- rbind(bp, mf, cc)
    res <- res[res$Adjp_BH < 0.05,]
    res <- res[order(res$Adjp_BH),]

    write.table(res, file = paste0("Signatures/results/data/GREAT/", arche, "_", analysis, ".tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
    return(res)
}
