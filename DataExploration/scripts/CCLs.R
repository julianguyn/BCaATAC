setwd("/Users/julianguyen/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(reshape2)
    library(ggpubr)
    library(grid)
    library(gridExtra)

})

source("source/DrugResponse/helper.R")
source("source/DataExploration/helper.R")
source("source/DataExploration/plots.R")
source("source/palettes.R")

###########################################################
# Load in data
###########################################################

# load in signature scores using sourced get_all_scores()
signature_scores <- get_all_scores()

# drug sensitivity processed from DrugResponse/ARCHE_CCLs.R
load("DrugResponse/results/data/sensitivity_data.RData")

###########################################################
# Plot distribution of ARCHE scores
###########################################################

# format dataframe to plot
toPlot <- t(signature_scores) |> as.data.frame()
toPlot <- data.frame(ARCHE = rep(paste0('ARCHE', 1:6), each = nrow(toPlot)),
                     Value = c(toPlot$ARCHE1, toPlot$ARCHE2, toPlot$ARCHE3,
                               toPlot$ARCHE4, toPlot$ARCHE5, toPlot$ARCHE6))
# plot distribution
plot_ARCHE_scores(toPlot)

###########################################################
# Get cell line overlap
###########################################################

all_cl <- unique(c(
    colnames(ubr1_sen),
    colnames(ubr2_sen),
    colnames(gray_sen),
    colnames(gcsi_sen),
    colnames(gdsc_sen),
    colnames(ctrp_sen),
    colnames(ccle_sen)
))

toPlot <- data.frame(
    UBR1 = ifelse(all_cl %in% colnames(ubr1_sen), 1, 0),
    UBR2 = ifelse(all_cl %in% colnames(ubr2_sen), 1, 0),
    GRAY = ifelse(all_cl %in% colnames(gray_sen), 1, 0),
    gCSI = ifelse(all_cl %in% colnames(gcsi_sen), 1, 0),
    GDSC2 = ifelse(all_cl %in% colnames(gdsc_sen), 1, 0),
    CTRP = ifelse(all_cl %in% colnames(ctrp_sen), 1, 0),
    CCLE = ifelse(all_cl %in% colnames(ccle_sen), 1, 0))
toPlot$sample <- all_cl
toPlot <- melt(toPlot)
toPlot$variable <- factor(toPlot$variable, levels = names(PSet_pal))

# plot overlap
plot_BCa_CCLs_overlap(toPlot)

###########################################################
# Get drug overlap
###########################################################

ubr1 <- drug_overlap(ubr1_sen, "UBR1")
ubr2 <- drug_overlap(ubr2_sen, "UBR2")
gray <- drug_overlap(gray_sen, "GRAY")
gcsi <- drug_overlap(gcsi_sen, "gCSI")
gdsc <- drug_overlap(gdsc_sen, "GDSC2")
ctrp <- drug_overlap(ctrp_sen, "CTRP")
ccle <- drug_overlap(ccle_sen, "CCLE")

###########################################################
# Plot BCa drugs overlapping in PSets
###########################################################

# compile results
toPlot <- rbind(ubr1, ubr2, gray, gcsi, gdsc, ctrp, ccle)
toPlot$Present <- factor(toPlot$Present, levels = c("Present", "Absent"))
toPlot$Drug <- factor(toPlot$Drug, levels = rev(bca_drugs))

# plot overlap
plot_BCa_drugs_overlap(toPlot)

###########################################################
# Assess drug response correlation across PSets
###########################################################

# create dataframe to hold correlation results
corr_res <- data.frame(matrix(nrow=0, ncol=3))
colnames(corr_res) <- c("PSet1", "PSet2", "PCC")

# plot for all drugs
corr_res <- plot_allDrugCorr(ubr1_sen, ubr2_sen, gray_sen, gcsi_sen, gdsc_sen, ctrp_sen, ccle_sen, "All", corr_res)

# same thing but subset for bca drugs
corr_res <- data.frame(matrix(nrow=0, ncol=3))
colnames(corr_res) <- c("PSet1", "PSet2", "PCC")
corr_res <- plot_allDrugCorr(
    ubr1_sen[rownames(ubr1_sen) %in% bca_drugs,],
    ubr2_sen[rownames(ubr2_sen) %in% bca_drugs,],
    gray_sen[rownames(gray_sen) %in% bca_drugs,],
    gcsi_sen[rownames(gcsi_sen) %in% bca_drugs,],
    gdsc_sen[rownames(gdsc_sen) %in% bca_drugs,],
    ctrp_sen[rownames(ctrp_sen) %in% bca_drugs,],
    ccle_sen[rownames(ccle_sen) %in% bca_drugs,],
    "BCa", corr_res)
