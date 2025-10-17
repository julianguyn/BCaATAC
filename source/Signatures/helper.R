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
