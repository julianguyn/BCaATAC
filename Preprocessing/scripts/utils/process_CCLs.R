#' Load in RNA-Seq counts matrix from other PSets
#' 
#' @param input string. Path to input PSet
#' @param output string. Path to save dataframe
#' 
get_pset_rna <- function(input, output) {
    pset <- readRDS(input) |> updateObject()
    pset <- summarizeMolecularProfiles(pset, mDataType = "Kallisto_0.46.1.rnaseq.counts")
    rna <- pset@assays@data$expr

    # keep only cell lines being used
    samples <- get_cells()
    rna <- rna[,colnames(rna) %in% samples$sample]
    df <- cbind(data.frame(Genes = rownames(rna), rna))
    colnames(df) <- c("Genes", colnames(rna))

    # save
    write.table(df, file = output,  quote = F, sep = "\t", col.names = T, row.names = F)
}
