#' Get drug sensitivity 
#'
#' Extract drug sensitivity dataframe from PSet
#' @param pset string. Path to PSet from which to extract sensitivity data from
#' @param load boolean. If TRUE, load in PSet RDS object instead of using readRDS()
#' @param update boolean. If TRUE, use updateObject() on PSet loaded from readRDS()
#' @param map boolean. If TRUE, use mapping() to standardize sample nanmes
#' @return Dataframe of PSet drug sensitivities.
#' 
get_drugsen <- function(pset, load = FALSE, update = TRUE, map = FALSE) {

    # get cells
    meta <- read.csv("metadata/lupien_metadata.csv")
    meta <- meta[meta$type == "cell_line", ]

    # load in PSet
    if (load == FALSE) {
        if (update == TRUE) {
            pset <- readRDS(pset) |> updateObject()
        } else {
            pset <- readRDS(pset)
        } 
    } else {
        load(pset)
        pset <- CTRP
    }

    # get sensitivity data
    if (map == FALSE) {
        sen <- PharmacoGx::summarizeSensitivityProfiles(
            pset,
            sensitivity.measure = "aac_recomputed",
            fill.missing = FALSE
        ) |> as.data.frame()
    } else {
        # if map = TRUE, use map_sen() 
        sen <- dcast(pset@treatmentResponse$profiles, treatmentid ~ sampleid, value.var = "aac_recomputed")
        rownames(sen) <- sen$treatmentid
        sen$treatmentid <- NULL
        sen <- map_sen(sen)
        sen <- sen/100
    }
    rm(pset)

    # subset to only cells from samples
    if ("MCF-10A" %in% colnames(sen)) {
        colnames(sen)[colnames(sen) == "MCF-10A"] <- "MCF10A"
    }
    sen <- sen[,which(colnames(sen) %in% meta$sampleid)]
    sen <- sen[,order(colnames(sen))]
    message(paste("Number of cells:", ncol(sen)))
    return(sen)
}

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
