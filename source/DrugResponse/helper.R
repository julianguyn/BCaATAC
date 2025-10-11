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

#' Get unique 49 BCa cell lines for analysis
#'
#' For duplicated cell lines, keep only those from Nergiz
#' @return Subsetted dataframe of 49 sample filenames.
#' 
get_all_scores <- function() {

    # get samples
    samples <- get_cells()

    # load in signature scores
    signature_scores <- read.table("Signatures/data/bca_sign.Zscore.txt")

    # remove duplicates and rename signatures & cells
    signature_scores <- signature_scores[,which(colnames(signature_scores) %in% samples$file)]
    rownames(signature_scores) <- paste0("ARCHE", 1:6)
    colnames(signature_scores) <- samples$sample[match(colnames(signature_scores), samples$file)] 
    signature_scores <- signature_scores[,order(colnames(signature_scores))]

    return(signature_scores)
}

#' Standardize PSet cell and drug naming
#'
#' Reference defined mapping to standardize PSet cell and drug names.
#' @param sen dataframe. Sensitivity dataframe obtained from get_drugsen()
#' @return Sen dataframe with standardized cell and drug names.
#' 

# define mappings
mapping_cells <- c("BT20" = "BT-20", 
                    "BT474" = "BT-474", 
                    "BT549" = "BT-549", 
                    "CAL120" = "CAL-120", 
                    "CAL148" = "CAL-148", 
                    "CAL51" = "CAL-51", 
                    "CAMA1" = "CAMA-1", 
                    "EFM19" = "EFM-19", 
                    "HDQP1" = "HDQ-P1", 
                    "HS578T" = "Hs 578T",
                    "JIMT1" = "JIMT-1", 
                    "KPL1" = "KPL-1", 
                    "LY2" = "MCF-7/LY2", 
                    "MCF7" = "MCF-7", 
                    "MDAMB157" = "MDA-MB-157",
                    "MDAMB175VII" = "MDA-MB-175-VII", 
                    "MDAMB231" = "MDA-MB-231", 
                    "MDAMB361" = "MDA-MB-361", 
                    "MDAMB436" = "MDA-MB-436", 
                    "MDAMB468" = "MDA-MB-468", 
                    "MX1" = "MX-1", 
                    "SKBR3" = "SK-BR-3", 
                    "SKBR5" = "SK-BR-5", 
                    "SKBR7" = "SK-BR-7", 
                    "SUM149" = "SUM149PT", 
                    "SUM159" = "SUM159PT", 
                    "T47D" = "T-47D",
                    "UACC812" = "UACC-812", 
                    "ZR751" = "ZR-75-1")
mapping_drugs <- c("carboplatin" = "Carboplatinum", 
                    "dipyridamole" = "Dipyridamole", 
                    "epirubicin" = "Epirubicin", 
                    "eribulin" = "Eribulin", 
                    "fluvastatin" = "Fluvastatin", 
                    "herceptin" = "Trastuzumab", 
                    "lapatinib" = "Lapatinib", 
                    "paclitaxel" = "Paclitaxel",
                    "OLA" = "Olaparib")

# function to standardize mapping of sensitivity dataframe
map_sen <- function(sen) {

    for (i in 1:length(colnames(sen))) {
        cell = colnames(sen)[i]
        if (cell %in% names(mapping_cells)) {colnames(sen)[i] <- unname(mapping_cells[cell])}
    }

    for (i in 1:length(rownames(sen))) {
        cell = rownames(sen)[i]
        if (cell %in% names(mapping_drugs)) {rownames(sen)[i] <- unname(mapping_drugs[cell])}
    }
    return(sen)
}

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

    # get 49 samples
    samples <- get_cells()

    # load in PSet
    if (load == FALSE) {
        if (update == TRUE) {
            pset <- readRDS(pset) |> updateObject()
        } else {
            pset <- readRDS(pset)
        } 
    } else {
        load(pset)
        pset = CTRP
    }

    # get sensitivity data
    if (map == FALSE) {
        sen <- PharmacoGx::summarizeSensitivityProfiles(pset, 
                                                      sensitivity.measure = "aac_recomputed",  
                                                      fill.missing = F) |> as.data.frame()
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
    sen <- sen[,which(colnames(sen) %in% samples$sample)]
    sen <- sen[,order(colnames(sen))]
    message(paste("Number of cells:", ncol(sen)))
    return(sen)
}

#' Get signature scores
#'
#' Subset signature scores for available cell lines for each PSet
#' @param sen dataframe. Drug sensitivity df from get_drugsen() for given PSet
#' @return Dataframe of PSet-specific signature scores.
#' 
get_scores <- function(sen) {
    to_keep <- samples[which(samples$sample %in% colnames(sen)),]
    sig <- signature_scores[,which(colnames(signature_scores) %in% to_keep$sample)]
    sig <- sig[,order(colnames(sig))]
    return(sig)
}

#' FDR correction and format dataframes from computeCI() and computePC()
#'
#' @param combinations dataframe. Intermediary from computeCI() and computePC().
#' @return A dataframe with FDR correction and additional formating.
#' 
format_combinations <- function(combinations, label) {
    
    # filtering and multiple test correction
    combinations <- combinations[complete.cases(combinations$pvalue),]
    combinations$FDR <- p.adjust(combinations$pvalue, method = "BH", n = length(combinations$pvalue))
    combinations$FDRsig <- ifelse(combinations$FDR < 0.05, TRUE, FALSE)

    # format dataframe for plotting (ordered by CI/PC already)
    combinations$rank <- 1:nrow(combinations)
    combinations$pairs <- paste0(combinations$signature, "_", combinations$drug)
    combinations$pset <- c(rep(label, nrow(combinations)))

    return(combinations)
}

#' Compute concordance index
#'
#' Computes concordance index between ARCHE scores and drug AAC.
#' @param signature_scores dataframe. ARCHE scores from get_scores().
#' @param sensitivity_data dataframe. Drug sensitivity from get_drugsen().
#' @param label string. PSet name
#' @return A dataframe of CI and metrics.
#' 
computeCI <- function(signature_scores, sensitivity_data, label) {

    # create data frame to hold results
    combinations <- as.data.frame(matrix(data = NA, nrow = nrow(signature_scores) * nrow(sensitivity_data), ncol = 7))
    colnames(combinations) <- c("signature", "drug", "ci", "pvalue", "se", "upper", "lower")
    combinations$signature <- rep(rownames(signature_scores), nrow(sensitivity_data))
    combinations$drug <- rep(rownames(sensitivity_data), each = nrow(signature_scores))

    # compute concordance index
    for (i in 1:nrow(combinations)){

        ci <- survcomp::concordance.index(as.numeric(sensitivity_data[combinations$drug[i],]), # drug AAC
                        surv.time = as.numeric(unlist(-signature_scores[combinations$signature[i],])), # ARCHE score
                        surv.event = rep(1,length(sensitivity_data)), 
                        outx = TRUE, method="noether", na.rm = TRUE)

        combinations$pvalue[i] <- ci$p.value
        combinations$ci[i] <- ci$c.index
        combinations$se[i] <- ci$se
        combinations$upper[i] <- ci$upper
        combinations$lower[i] <- ci$lower
    }

    # format dataframe and FDR correction
    combinations <- combinations[order(combinations$ci),]
    combinations <- format_combinations(combinations, label)
    
    return(combinations)
}

#' Compute Pearson's correlation
#'
#' Computes PC between ARCHE scores and drug AAC.
#' @param signature_scores dataframe. ARCHE scores from get_scores().
#' @param sensitivity_data dataframe. Drug sensitivity from get_drugsen().
#' @param label string. PSet name
#' @return A dataframe of CI and metrics.
#' 
computePC <- function(signature_scores, sensitivity_data, label) {
    
    # create data frame to hold results
    combinations <- as.data.frame(matrix(data = NA, nrow = nrow(signature_scores) * nrow(sensitivity_data), ncol = 4))
    colnames(combinations) <- c("signature", "drug", "pc", "pvalue")
    combinations$signature <- rep(rownames(signature_scores), nrow(sensitivity_data))
    combinations$drug <- rep(rownames(sensitivity_data), each = nrow(signature_scores))

    # compute concordance index
    for (i in 1:nrow(combinations)){
        AAC <- as.numeric(sensitivity_data[combinations$drug[i],])
        if (length(AAC[!is.na(AAC)]) > 2) {
            pc <- cor.test(AAC, # drug AAC
                    as.numeric(unlist(signature_scores[combinations$signature[i],])), # ARCHE score
                    method = 'pearson', alternative = 'two.sided')

            combinations$pvalue[i] <- pc$p.value
            combinations$pc[i] <- pc$estimate
        }
    }

    # format dataframe and FDR correction
    combinations <- combinations[order(combinations$pc),]
    combinations <- format_combinations(combinations, label)
    
    return(combinations)
}

#' Extract AAC for drug of interest
#'
#' Helper function for plot_CCLs().
#' @param pset dataframe. PSet drug sensitivity from get_drugsen()
#' @param drugs string. Drug of interest.
#' @param df dataframe. Dataframe to bind results to.
#' @param label string. PSet name.
#' @return df with drug AAC binded (if present in PSet).
#' 
get_doiAAC <- function(pset, drug, df, label) {
    if (drug %in% rownames(pset)) {
        # save drug response
        res <- melt(pset[rownames(pset) == drug,]) |> suppressMessages()
        colnames(res) <- c("Sample", "AAC")
        res$PSet <- label
        # merge results with compiled dataframe
        df <- rbind(df, res)
        return(df)
    } else {return(df)}
}

#' Save significant concordance index results
#'
#' Append columns to treatment response matrices with PDX ARCHE signature scores.
#' @param df string. Vector of PDX model names standardize against mapping.
#' @return A vector of (string) standardized PDX model names.
#' 
saveSig <- function(sig_com, combinations, label) {

    # save signatures associated with drug sensitivity:
    sen <- combinations[which(combinations$ci > 0.6 & combinations$FDR < 0.05),]
    sig_com <- rbind(sig_com, data.frame(pset = c(rep(label, nrow(sen))), signature = sen$signature, drug = sen$drug, ci = sen$ci, pvalue = sen$pvalue, FDR = sen$FDR))
    #write.csv(sen[order(sen$FDR, -sen$ci),], file = paste0("DrugResponse/results/tables/indiv_PSet_CI/", label, "_sen.csv"), row.names = F)

    # save signatures associated with drug resistance:
    res <- combinations[which(combinations$ci < 0.4 & combinations$FDR < 0.05),]
    sig_com <- rbind(sig_com, data.frame(pset = c(rep(label, nrow(res))), signature = res$signature, drug = res$drug, ci = res$ci, pvalue = res$pvalue, FDR = res$FDR))
    #write.csv(res[order(res$FDR, res$ci),], file = paste0("DrugResponse/results/tables/indiv_PSet_CI/", label, "_res.csv"), row.names = F)

    return(sig_com)
}