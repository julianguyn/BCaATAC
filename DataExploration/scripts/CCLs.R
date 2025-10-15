setwd("/Users/julianguyen/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
})

source("source/DrugResponse/helper.R")
source("source/DataExploration/helper.R")
source("source/DataExploration/plots.R")

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
