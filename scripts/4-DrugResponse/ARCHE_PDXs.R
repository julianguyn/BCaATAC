# load libraries
suppressPackageStartupMessages({
    library(readxl)
    library(ggplot2)
    library(RColorBrewer)
    library(dplyr)
    library(ROCR)
    library(ggpubr)
    library(ggh4x)
    library(grid)
    library(gridExtra)
})

source("utils/plots/drug_response_pdx.R")
source("utils/plots/drug_response_pdx_indivplots.R")
source("utils/palettes.R")
source("utils/get_data.R")
source("utils/mappings.R")

###########################################################
# Load in data
###########################################################

# pdx drug response data
xeva1 <- get_xeva("data/rawdata/pdx/DrugResponse_PDX.xlsx", xlsx = TRUE)
xeva2 <- get_xeva_july()
xeva3 <- get_xeva("data/rawdata/pdx/PDX_Response_Sept2025.csv")

# get all ARCHE scores
scores <- get_arche_pdx()

###########################################################
# Assign signature scores
###########################################################

# helper function
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

# get ARCHE scores per xeva collection
xeva1_sig <- get_ARCHE(xeva1)
xeva2_sig <- get_ARCHE(xeva2)
xeva3_sig <- get_ARCHE(xeva3)

###########################################################
# Assess ARCHE drug response associations 
###########################################################

x1 <- assess_ARCHE_PDX(xeva1_sig, "Xeva1") |> suppressWarnings()
x2 <- assess_ARCHE_PDX(xeva2_sig, "Xeva2") |> suppressWarnings()
x3 <- assess_ARCHE_PDX(xeva3_sig, "Xeva3")

##todo: format for x3: no AUC! mre, br, and bar
##todo: why does it keep printing NULL out