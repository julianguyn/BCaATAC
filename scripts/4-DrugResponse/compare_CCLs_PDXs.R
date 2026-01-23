
###########################################################
# Load in data
###########################################################

# load in PDX results
path <- "data/results/data/4-DrugResponse/PDX/rmNergizKomal-fullXeva/"
files <- list.files(path, full.names = TRUE)

df <- data.frame(matrix(nrow=0, ncol=0))
for (file in files) {
    res <- read.csv(file)
    df <- rbind(df, res)
}

# load in cell line results
path <- "data/results/data/4-DrugResponse/CCLs/"
files <- list.files(path, full.names = TRUE)
files <- files[grep("ClassA_Biomarkers_", files)]

cl <- data.frame(matrix(nrow=0, ncol=0))
for (file in files) {
    res <- read.csv(file)
    res$ARCHE_label <- sub(".*Biomarkers_", "", sub("\\.csv", "", file))
    cl <- rbind(cl, res)
}

###########################################################
# Explore results
###########################################################

# set thresholds
pc <- 0.4
pval <- 0.1

# get sig pdx drugs
sig_bar <- df[df$PC.BAR_median > pc & df$pval.BAR_median < pval,]
pdx_drugs <- unique(sig_bar$drug)

# get sig ccl drugs
ccl_drugs <- unique(cl$drug)

###########################################################
# Define matches
###########################################################

top1i_pdx <- c("TOPOTECAN", "TOPOTECAN-LOW", "SACITUZUMAB-GOVITECAN", "DERUXTECAN")
top1i_ccl <- c("Topotecan", "SN-38")
top2i_ccl <- c("Etoposide", "Teniposide")

pactx_pdx <- pdx_drugs[grep("PACLITAXEL", pdx_drugs)]
pactx_ccl <- "Paclitaxel"

trast_pdx <- pdx_drugs[grep("TRASTUZUMAB", pdx_drugs)]
trast_ccl <- "Trastuzumab"

bcl_2_pdx <- pdx_drugs[grep("ABT-263", pdx_drugs)]
bcl_2_ccl <- c("Navitoclax", ccl_drugs[grep("navitoclax", ccl_drugs)])

parpi_pdx <- c("AZD-5305")
parpi_ccl <- c("Olaparib", "Veliparib")

statn_pdx <- pdx_drugs[grep("FLUVASTATIN", pdx_drugs)]
statn_ccl <- c("Fluvastatin", "Simvastatin", "Lovastatin")

dipyr_pdx <- c("DIPYRIDAMOLE")
dipyr_ccl <- c("Dipyridamole")

sig_bar[sig_bar$drug %in% top1i_pdx,]

###########################################################
# Search matches
###########################################################

# helper function
search_match <- function(pdx_drug, ccl_drug) {

    pdx_df <- sig_bar[sig_bar$drug %in% pdx_drug,]
    arches <- unique(pdx_df$ARCHE)

    ccl_df <- cl[cl$drug %in% ccl_drug,]
    match <- FALSE
    for (arche in arches) {
        if (arche %in% ccl_df$signature) {
            print(paste("Match:", arche))
            print(ccl_df[ccl_df$signature == arche,c(1:3,8,11:13)])
            match <- TRUE
        }
    }
    if (match == FALSE) print(paste("\nNo match found for", pdx_drug))
}

#search_match(top1i_pdx, top1i_ccl)
search_match(top1i_pdx, top2i_ccl) #ARCHE4,6
search_match(pactx_pdx, pactx_ccl) #ARCHE2,5,6
#search_match(trast_pdx, trast_ccl)
search_match(bcl_2_pdx, bcl_2_ccl) #ARCHE2,5
search_match(parpi_pdx, parpi_ccl) #ARCHE1
search_match(statn_pdx, statn_ccl) #ARCHE2,3,4,5,6
#search_match(dipyr_pdx, dipyr_ccl)