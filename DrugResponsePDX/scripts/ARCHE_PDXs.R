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

source("source/DrugResponsePDX/helper.R")
source("source/DrugResponsePDX/plots.R")
source("source/DrugResponsePDX/indiv_plots.R")
source("source/palettes.R")

###########################################################
# Load in drug response data
###########################################################

xeva1 <- get_xeva("DrugResponsePDX/data/drugresponse/DrugResponse_PDX.xlsx", xlsx = TRUE)
xeva2 <- get_xeva_july()
xeva3 <- get_xeva("DrugResponsePDX/data/Sept2025/PDX_Response_Sept2025.csv")

###########################################################
# Assign signature scores
###########################################################

# get all ARCHE scores
scores <- get_pdx_ARCHE()

# quick check for mismatched samples
check_sample_overlap(rownames(scores), xeva1$patient.id, "ARCHE Scores", "Xeva1")
check_sample_overlap(rownames(scores), xeva2$patient.id, "ARCHE Scores", "Xeva2")
check_sample_overlap(rownames(scores), xeva3$patient.id, "ARCHE Scores", "Xeva3")

# get ARCHE scores per xeva collection
xeva1_sig <- get_ARCHE(xeva1)
xeva2_sig <- get_ARCHE(xeva2)
xeva3_sig <- get_ARCHE(xeva3)

###########################################################
# Assess ARCHE drug response associations 
###########################################################

x2 <- assess_ARCHE_PDX(xeva2_sig, "Xeva2")