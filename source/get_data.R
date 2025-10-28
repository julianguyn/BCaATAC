#' Get unique 49 BCa cell lines for analysis
#'
#' For duplicated cell lines, keep only those from Nergiz
#' @return Subsetted dataframe of 49 sample filenames.
#' 
get_cells <- function() {

    # load in BCa cell lines
    samples <- read.csv("DrugResponse/data/cl_samples.csv")

    # get duplicates
    dup <- samples$sample[duplicated(samples$sample)]
    tmp <- samples[samples$sample %in% dup,]

    # keep Nergiz dups
    samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),])

    return(samples)
}

#' Get CCL gene counts matrices
#' 
#' @param pset string. Name of PSet to load
#' 
get_pset_rna <- function(pset) {

    pset_options <- c("UBR2", "GRAY", "gCSI", "CCLE")
    if (!pset %in% pset_options) {
        message("Option not valid. Pick from: UBR2, GRAY, gCSI, CCLE")
        return(NA)
    }

    filename <- paste0("Preprocessing/procdata/CCLs/", pset, "_RNA.tsv")
    df <- fread(filename, data.table = FALSE)
    df$Genes <- gsub("\\..*", "", df$Genes)
    df <- df %>%
        dplyr::group_by(Genes) %>%
        dplyr::summarise(across(everything(), mean)) %>%
        as.data.frame()
    rownames(df) <- df$Genes
    df$Genes <- NULL
    return(df)
}
