# load libraries
suppressPackageStartupMessages({
    library(readxl)
    library(ggplot2)
    library(RColorBrewer)
    library(dplyr)
    library(ROCR)
    library(ggpubr)
})

source("source/DrugResponsePDX/helper.R")

###########################################################
# Load in data
###########################################################

# read in drug response data
old_response <- as.data.frame(read_excel("DrugResponsePDX/data/drugresponse/DrugResponse_PDX.xlsx", sheet = 1))
new_xeva_auc <- read.csv("DrugResponsePDX/data/Jul232025-Xeva/auc_reps.csv")
new_xeva_mRECIST <- read.csv("DrugResponsePDX/data/Jul232025-Xeva/mRECIST_reps.csv")
bar <- as.data.frame(read_excel("DrugResponsePDX/data/Jul232025-Xeva/PDX_BR_BAR_grouped_Julia_8205_TAX_selectedmodels.xlsx", sheet = 1))
bar$patient.id <- bar$PDX_ID

# read in signature scores
scores <- as.data.frame(t(read.table("DrugResponsePDX/data/chromvar/bca_sign.Zscore.txt")))
colnames(scores) <- paste0("Signature", 1:6)


###########################################################
# Format sample names
###########################################################

# using sourced map_pdx() to standardize PDX names
rownames(scores) <- map_pdx(rownames(scores))
old_response$patient.id <- map_pdx(old_response$patient.id)
new_xeva_auc$patient.id <- map_pdx(new_xeva_auc$patient.id)
new_xeva_mRECIST$patient.id <- map_pdx(new_xeva_mRECIST$patient.id)
bar$patient.id <- map_pdx(bar$patient.id)


###########################################################
# Assign signature scores
###########################################################

# using sourced get_ARCHE() to extract PDX ARCHE scores
old_response <- get_ARCHE(old_response)
new_xeva_auc <- get_ARCHE(new_xeva_auc)
new_xeva_mRECIST <- get_ARCHE(new_xeva_mRECIST)
bar <- get_ARCHE(bar)

###########################################################
# Assess ARCHE drug response associations (via mRECIST)
###########################################################

# using sourced assess_ARCHE_mRECIST()
assess_ARCHE_mRECIST(old_response, 'AZD-5305')
assess_ARCHE_mRECIST(old_response, 'AZD-8205')
assess_ARCHE_mRECIST(old_response, 'BMN-673')
assess_ARCHE_mRECIST(old_response, 'CARBOPLATIN')
assess_ARCHE_mRECIST(old_response, 'CFI-400945')
assess_ARCHE_mRECIST(old_response, 'CFI-402257')
assess_ARCHE_mRECIST(old_response, 'ERIBULIN')
assess_ARCHE_mRECIST(old_response, 'TAXOL')

assess_ARCHE_mRECIST(new_xeva_mRECIST, 'AZD-5305')
assess_ARCHE_mRECIST(new_xeva_mRECIST, 'AZD-8205')
assess_ARCHE_mRECIST(new_xeva_mRECIST, 'CDX-011')
assess_ARCHE_mRECIST(new_xeva_mRECIST, 'CFI-400945')
assess_ARCHE_mRECIST(new_xeva_mRECIST, 'CX-5461')
assess_ARCHE_mRECIST(new_xeva_mRECIST, 'FLUVASTATIN')

###########################################################
# Assess ARCHE drug response associations (via AUC)
###########################################################

assess_ARCHE_AUC(old_response, 'AZD-5305')
assess_ARCHE_AUC(old_response, 'AZD-8205')
assess_ARCHE_AUC(old_response, 'BMN-673')
assess_ARCHE_AUC(old_response, 'CARBOPLATIN')
assess_ARCHE_AUC(old_response, 'CFI-400945')
assess_ARCHE_AUC(old_response, 'CFI-402257') 
assess_ARCHE_AUC(old_response, 'ERIBULIN')
assess_ARCHE_AUC(old_response, 'TAXOL')

assess_ARCHE_AUC(new_xeva_auc, 'AZD-5305')
assess_ARCHE_AUC(new_xeva_auc, 'AZD-8205')
assess_ARCHE_AUC(new_xeva_auc, 'CDX-011')
assess_ARCHE_AUC(new_xeva_auc, 'CFI-400945')
assess_ARCHE_AUC(new_xeva_auc, 'CX-5461')
assess_ARCHE_AUC(new_xeva_auc, 'FLUVASTATIN')

###########################################################
# Assess ARCHE drug response associations (via BAR)
###########################################################

assess_ARCHE_BAR(bar, 'AZD-8205')
assess_ARCHE_BAR(bar, 'PACLITAXEL')