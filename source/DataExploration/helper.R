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