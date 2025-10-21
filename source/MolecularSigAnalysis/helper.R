#' Run DEG
#' 
#' @param counts dataframe. Counts matrix
#' @param meta dataframe. Metadata table
#' @param arche string. ARCHE to run by
#' @param subtype string. Subtype to run by 
#' 
run_DEG <- function(counts, meta, arche = NA, subtype = NA) {

  # create colData
  colData <- data.frame(
    SampleID = meta$Sample.Name,
    ARCHE = ifelse(meta$Signature == arche, arche, "Other"),
    Subtype = ifelse(meta$Subtype == subtype, subtype, "Other")
  )
  
  if (!is.na(arche)) {  # if running ARCHEx vs all
    colData$ARCHE <- factor(colData$ARCHE, levels = c("Other", arche))
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design= ~ ARCHE) |>
                                suppressMessages()
    label <- paste0("ARCHE_", arche, "_vs_Other")
  } else {            # if running Subtypex vs all
    colData$Subtype <- factor(colData$Subtype, levels = c("Other", subtype))
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design= ~ Subtype) |>
                                suppressMessages()
    label <- paste0("Subtype_", subtype, "_vs_Other")
  }

  # DEG
  dds <- DESeq(dds) |> suppressMessages()
  res <- results(dds, name=label)

  return(res)
}
