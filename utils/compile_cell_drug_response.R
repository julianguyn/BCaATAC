# -----
# Load in drug match lists
source("utils/pdx_drug_matches.R")

#' Helper function to compile cell drug response
#' associations with same drug or drugs of similar 
#' MoA as list provided
#' 
#' Usage: validate PDX associations
#' 
#' @param pair string. In format of "ARCHEx_DRUG"
#' @param df data.frame. Cell drug response associations
#' with columns 'signature', 'drug', 'pc', 'pset', 'FDR', and 'label'
#' @param drug_list list. List from utils/pdx_drug_matches.R
#' 
search_drug_response_cells <- function(pair, df, drug_list = drug_list1) {
    arche <- sub("_.*", "", pair)
    drug_to_search <- sub(".*_", "", pair)
    matched_drugs <- drug_list[[drug_to_search]]

    matches <- df[
        which(df$signature == arche & df$drug %in% matched_drugs & df$FDR < 0.1 & abs(df$pc) > 0.39),
        c("drug", "pset", "pc", "FDR", "Label")
    ]
    rownames(matches) <- NULL

    message("----------------------------------------")
    message(paste("Search results for", pair))

    for (drug in matched_drugs) {
        message(paste("\n", arche, "and", drug))
        toPrint <- matches[matches$drug == drug,]
        toPrint <- toPrint[order(abs(toPrint$pc), decreasing = TRUE),]
        print(toPrint)
    }
    message("\n")
}

#  search_drug_response_cells("ARCHE5_PACLITAXEL", cell_toPlot)

#' Helper function to compile drug response
#' associations for all drugs from lists
#' 
#' Usage: quick search
#' 
#' @param arche string. ARCHE to search
#' @param df data.frame. Cell drug response associations
#' with columns 'signature', 'drug', 'pc', 'pset', 'FDR', and 'label'
#' 
validate_drug_response_cells <- function(arche, df) {
    for (drug in names(drug_list2)) {
        pair <- paste0(arche, "_", drug)
        search_drug_response_cells(pair, df, drug_list2)
    }
}

#  validate_drug_response_cells("ARCHE5", df)

load(paste0("data/results/data/Misc/cell_pdx_tcga_ARCHE_drug_response_testing_all_scoring.RData"))

### --- SANDBOX
# CFI-400945
# PACLITAXEL-15DAILY
# ERIBULIN-LOW
# TRASTUZUMAB-DERUXTECAN

search_drug_response_cells("ARCHE6_EVEROLIMUS", cell_toPlot)
validate_drug_response_cells("ARCHE5", cell_toPlot)
