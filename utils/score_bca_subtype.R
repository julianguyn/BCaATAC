#' Compute BCa subtypes using genefu
#' 
#' @param rna dataframe. RNA-Seq counts matrix
#' @param meta dataframe. Gene metadata
#' @param model string. Subtyping classification model for genefu
#' Options: "scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS", or "claudinLow"
#' @return dataframe of subtype scores per sample
score_bca_subtype <- function(rna, meta, model) {

    # format matrix
    rna <- t(rna) |> as.data.frame()

    # set do.mapping param
    do.mapping <- ifelse(model == "pam50", FALSE, TRUE)

    if (do.mapping == TRUE) {   # if not PAM50
        subset_meta <- meta[-which(is.na(meta$EntrezGene.ID)),]
        rna <- rna[,match(subset_meta$Ensembl, colnames(rna))]
        colnames(rna) <- subset_meta$EntrezGene.ID[match(subset_meta$Ensembl, colnames(rna))]

        # average across duplicates
        unique_genes <- unique(colnames(rna))
        rna <- as.data.frame(sapply(unique_genes, function(x) {
            rowMeans(rna[, names(rna) == x, drop = FALSE])
        }))
        meta <- subset_meta[!duplicated(subset_meta$EntrezGene.ID),]
        dimnames(meta)[[1]] <- meta$EntrezGene.ID
    } else {                    # if PAM50
        # match and get gene names
        colnames(rna) <- meta$Gene.Symbol[match(colnames(rna), meta$Ensembl)]
    }

    # get subtype predictions for PAM50
    SubtypePredictions <- molecular.subtyping(sbt.model = model, 
                                              data = rna,
                                              annot = meta, 
                                              do.mapping = do.mapping)

    # reutrn subtypes and subtyping probability
    res_df <- cbind(Subtype = as.character(SubtypePredictions$subtype), 
                    as.data.frame(SubtypePredictions$subtype.proba))
    return(res_df)
}

# TODO: fix for other psets, try other methods