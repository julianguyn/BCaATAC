#' List of BCa-relevant drugs
#'
bca_drugs <- c("5-Fluorouracil",
               "Abemaciclib",
               "Alpelisib",
               "Anastrozole",
               "AZD-5305",
               "AZD-8205",
               "BMN-673",
               "Capecitabine",
               "Carboplatinum",
               "CDX-011",
               "CFI-400945",
               "CFI-402257",
               "Cisplatin",
               "CX-5461",
               "Cyclophosphamide",
               "Docetaxel",
               "Doxorubicin",
               "Elacestrant",
               "Epirubicin",
               "Eribulin",
               "Erlotinib",
               "Everolimus",
               "Exemestane",
               "Fluvastatin",
               "Fulvestrant",
               "Gemcitabine",
               "Inavolisib",
               "Ipatasertib",
               "Lapatinib",
               "Letrozole",
               "Neratinib",
               "Olaparib",
               "Paclitaxel",
               "Palbociclib",
               "Pembrolizumab",
               "Pertuzumab",
               "Ribociclib",
               "Selumetinib",
               "Tamoxifen",
               "Trastuzumab",
               "Tucatinib",
               "UNC0642")

#' Search for BCa drugs in PSets
#'
#' Searches for BCa drugs in pset drug sensitivity data.
#' @param sen dataframe. Drug sensitivity dataframe
#' @param label string. PSet label
#' @return A dataframe of meta analysis results.
#' 
drug_overlap <- function(sen, label) {
    tested_drugs <- rownames(sen)
    df <- data.frame(Drug = bca_drugs, Present = "Absent")

    for (i in 1:length(bca_drugs)) {
        drug = bca_drugs[i]
        if (drug %in% tested_drugs) {
            df$Present[i] <- "Present"
        }
    }
    overlap <- tested_drugs[tested_drugs %in% bca_drugs]
    print(paste("Number of overlapping drugs:", length(overlap)))
    df$PSet <- label
    return(df)
}

#' Creates dataframe for drug correlation and plotting
#'
#' @param pset1 dataframe. Drug sensitivity dataframe1
#' @param pset2 dataframe. Drug sensitivity dataframe2
#' @return dataframe object for correlation and plotting.
#' 
format_drug_pset <- function(pset1, pset2) {

    # get common samples
    common_cell <- intersect(colnames(pset1), colnames(pset2))
    common_drug <- intersect(rownames(pset1), rownames(pset2))
    pset1 <- pset1[match(common_drug, rownames(pset1)),match(common_cell, colnames(pset1))]
    pset2 <- pset2[match(common_drug, rownames(pset2)),match(common_cell, colnames(pset2))]

    # set labels and melt
    pset1$drug <- rownames(pset1)
    pset2$drug <- rownames(pset2)
    pset1 <- melt(pset1)
    pset2 <- melt(pset2)
    pset1$pairs <- paste0(pset1$variable, "_", pset1$drug)
    pset2$pairs <- paste0(pset2$variable, "_", pset2$drug)

    # # scatter plot of drug response difference for CTRP and GDSC 
    toPlot <- data.frame(
        pair = pset1$pairs, 
        pset1 = pset1$value, 
        pset2 = pset2$value[match(pset1$pairs, pset2$pairs)]
    )
    
    return(toPlot)
}