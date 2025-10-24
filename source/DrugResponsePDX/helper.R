#' Standardize PDX Naming
#'
#' Reference defined mapping to standardize PDX names.
#' @param models string. Vector of PDX model names standardize against mapping.
#' @return A vector of (string) standardized PDX model names.
#' 

# define mapping
mapping <- c("X108099P1" = "108099", 
             "108099P1" = "108099", 
             "X16720" = "16720",
             "X48602" = "48602",
             "X53782" = "53782",
             "BPTO95" = "BPTO.95",
             "INS_B014" = "INSB014",
             "INS_B019" = "INSB019",
             "REFS032" = "REF032",
             "REFS034P1" = "REF034", 
             "REFS038" = "REF038",
             "REFB047" = "REF047",
             "REF047_S26_L002" = "REF047",
             "REFS036P1A" = "REF036",
             "REF_S_038" = "REF038"
)

map_pdx <- function(models) {

    for (i in 1:length(models)) {
        model = models[i]
        if (model %in% names(mapping)) {
            models[i] <- unname(mapping[model])
        }
    }
    return(models)
}

#' Get PDX ARCHE scores
#'
#' @return ARCHE scores fataframe with standardized PDX IDs.
#' 
get_pdx_ARCHE <- function() {
    scores <- as.data.frame(t(read.table("DrugResponsePDX/data/chromvar/bca_sign.Zscore.txt")))
    colnames(scores) <- paste0("ARCHE", 1:6)
    rownames(scores) <- map_pdx(rownames(scores))
    return(scores)
}

#' Load in XEVA drug response spreadsheets
#' 
#' @param filepath string. 
#' @param xlsx boolean. 
#' 
get_xeva <- function(filepath, xlsx = FALSE) {

    # read in file
    if (xlsx == TRUE) {
        xeva <- as.data.frame(read_excel(filepath, sheet = 1))
    } else {
        xeva <- read.csv(filepath)
    }

    # create patient.id if needed and standardize
    if ("patient.id" %in% colnames(xeva)) {
        xeva$patient.id <- map_pdx(xeva$patient.id)
    } else if ("PDX_ID" %in% colnames(xeva)) {
        xeva$patient.id <- map_pdx(xeva$PDX_ID)
    } else {
        message("Check sampleID column name")
        return(NA)
    }

    # create model_group
    xeva$model_group <- paste(xeva$patient.id, xeva$drug, sep = "_")

    return(xeva)
}

#' Compile results from July Xeva data access
#' 
get_xeva_july <- function() {

    # load individiual files
    auc <- read.csv("DrugResponsePDX/data/Jul232025-Xeva/auc_reps.csv")
    mre <- read.csv("DrugResponsePDX/data/Jul232025-Xeva/mRECIST_reps.csv")
    bar <- as.data.frame(read_excel("DrugResponsePDX/data/Jul232025-Xeva/PDX_BR_BAR_grouped_Julia_8205_TAX_selectedmodels.xlsx", sheet = 1))

    # map sample names
    auc$patient.id <- map_pdx(auc$patient.id)
    mre$patient.id <- map_pdx(mre$patient.id)
    bar$patient.id <- map_pdx(bar$PDX_ID)

    # create model_group variable
    auc$model_group <- paste(auc$patient.id, auc$drug, sep = "_")
    mre$model_group <- paste(mre$patient.id, mre$drug, sep = "_")
    bar$model_group <- paste(bar$patient.id, bar$drug, sep = "_")

    # create new dataframe
    pairs <- unique(c(auc$model_group, mre$model_group, bar$model_group))
    xeva2 <- data.frame(
        model_group = pairs,
        patient.id =  gsub("_.*", "", pairs),
        drug = gsub(".*_", "", pairs)
    )
    xeva2$AUC <- auc$AUC[match(xeva2$model_group, auc$model_group)]
    xeva2$mRECIST <- mre$mRECIST[match(xeva2$model_group, mre$model_group)]
    xeva2$BR_median <- bar$BR_median[match(xeva2$model_group, bar$model_group)]
    xeva2$BAR_median <- bar$BAR_median[match(xeva2$model_group, bar$model_group)]

    # bar has same mRECIST as above but new mRECIST for paclitaxol (not in above)
    xeva2$mRECIST[xeva2$drug == "PACLITAXEL"] <- bar$mRECIST[match(xeva2$model_group[xeva2$drug == "PACLITAXEL"], bar$model_group)]
    return(xeva2)
}

#' Check overlap of sample IDs and report missing
#' 
#' Helper function for finding mismatched sample IDs lol
#' 
#' @param list1 string. First list of model IDs
#' @param list2 string. Second list of model IDs
#' @param label1 string. Name of first list
#' @param label2 string. Name of second list
#' 
check_sample_overlap <- function(list1, list2, label1, label2) {
    list1 <- unique(list1)
    list2 <- unique(list2)
    message(paste("\nNumber of samples from", label1, "in", label2))
    print(table(list1 %in% list2))
    message(paste("\nSamples from", label1, "NOT in", label2))
    missing <- list1[-which(list1 %in% list2)]
    print(missing[order(missing)])
    message(paste("\nNumber of samples from", label2, "in", label1))
    print(table(list2 %in% list1))
    message(paste("\nSamples from", label2, "NOT in", label1))
    missing <- list2[-which(list2 %in% list1)]
    print(missing[order(missing)])
}

#' Extract PDX ARCHE signature scores
#'
#' Append columns to treatment response matrices with PDX ARCHE signature scores.
#' @param df string. Vector of PDX model names standardize against mapping.
#' @return A vector of (string) standardized PDX model names.
#' 
get_ARCHE <- function(df) {
    df <- df[df$patient.id %in% rownames(scores),]
    df$ARCHE6 <- df$ARCHE5 <- df$ARCHE4 <- df$ARCHE3 <- df$ARCHE2 <- df$ARCHE1 <- NA
    for (i in 1:nrow(df)) {
        sample = df$patient.id[i]
        df$ARCHE1[i] <- scores[rownames(scores) == sample,]$ARCHE1
        df$ARCHE2[i] <- scores[rownames(scores) == sample,]$ARCHE2
        df$ARCHE3[i] <- scores[rownames(scores) == sample,]$ARCHE3
        df$ARCHE4[i] <- scores[rownames(scores) == sample,]$ARCHE4
        df$ARCHE5[i] <- scores[rownames(scores) == sample,]$ARCHE5
        df$ARCHE6[i] <- scores[rownames(scores) == sample,]$ARCHE6
        
    }
    return(df)
}
