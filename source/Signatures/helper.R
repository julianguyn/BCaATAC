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

#' Annotate peaks in ARCHE
#' 
#' @param gr GRanges object. GRanges of ARCHE peaks
#' @param arche string. ARCHE label
#' @return anno annotatePeak result
#' 
annotateARCHE <- function(gr, arche) {

    anno <- annotatePeak(gr, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
    anno$ARCHE <- arche
    return(anno)
}
