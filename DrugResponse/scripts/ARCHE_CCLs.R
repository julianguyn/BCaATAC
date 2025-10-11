setwd("/Users/julianguyen/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(survcomp)
    library(wesanderson)
    library(ggplot2)
    library(ggh4x)
    library(reshape2)
})

source("source/DrugResponse/helper.R")
source("source/DrugResponse/plots.R")
source("source/palettes.R")

###########################################################
# Load in data
###########################################################

# load in BCa cell lines 
samples <- get_cells()

# load in signature scores 
signature_scores <- get_all_scores()

###########################################################
# Load and save PSet sensitivity data
###########################################################

# get drug sensitivites
ubr1_sen <- get_drugsen("DrugResponse/data/PSet_UHNBreast.rds")                             #43 CCLs, 7 dups
ubr2_sen <- get_drugsen("DrugResponse/data/PharmacoSet.RDS", update = FALSE, map = TRUE)    #42 CCLs
gray_sen <- get_drugsen("DrugResponse/data/PSet_GRAY2017.rds")                              #42 CCLs
gcsi_sen <- get_drugsen("DrugResponse/data/gCSI.rds")                                       #25 CCLs, 7 dups
gdsc_sen <- get_drugsen("DrugResponse/data/GDSC2-8.2.rds")                                  #34 CCLs, 7 dups
ctrp_sen <- get_drugsen("DrugResponse/data/CTRP.rds", load = TRUE)                          #32 CCLs, 7 dups
ccle_sen <- get_drugsen("DrugResponse/data/CCLE.rds")                                       #23 CCLs, 7 dups

# save drug sensitivities
save(ubr1_sen, ubr2_sen, gray_sen, gcsi_sen, gdsc_sen, ctrp_sen, ccle_sen, 
     file = "DrugResponse/results/data/sensitivity_data.RData")

###########################################################
# Subset signature associations
###########################################################

# keep only CCLs with drug response
ubr1_sig <- get_scores(ubr1_sen)
ubr2_sig <- get_scores(ubr2_sen)
gray_sig <- get_scores(gray_sen)
gcsi_sig <- get_scores(gcsi_sen)
gdsc_sig <- get_scores(gdsc_sen)
ctrp_sig <- get_scores(ctrp_sen)
ccle_sig <- get_scores(ccle_sen)

###########################################################
# Compute CI and PC of ARCHE-drug associations
###########################################################

# PC
ubr1_PC <- computePC(ubr1_sig, ubr1_sen, "UBR1")
ubr2_PC <- computePC(ubr2_sig, ubr2_sen, "UBR2")
gray_PC <- computePC(gray_sig, gray_sen, "GRAY")
gcsi_PC <- computePC(gcsi_sig, gcsi_sen, "gCSI")
gdsc_PC <- computePC(gdsc_sig, gdsc_sen, "GDSC2")
ctrp_PC <- computePC(ctrp_sig, ctrp_sen, "CTRP")
ccle_PC <- computePC(ccle_sig, ccle_sen, "CCLE")

# CI
ubr1_CI <- computeCI(ubr1_sig, ubr1_sen, "UBR1")
ubr2_CI <- computeCI(ubr2_sig, ubr2_sen, "UBR2")
gray_CI <- computeCI(gray_sig, gray_sen, "GRAY")
gcsi_CI <- computeCI(gcsi_sig, gcsi_sen, "gCSI")
gdsc_CI <- computeCI(gdsc_sig, gdsc_sen, "GDSC2")
ctrp_CI <- computeCI(ctrp_sig, ctrp_sen, "CTRP")
ccle_CI <- computeCI(ccle_sig, ccle_sen, "CCLE")

# compile results
PC_res <- rbind(ubr1_PC, ubr2_PC, gray_PC, gcsi_PC, gdsc_PC, ctrp_PC, ccle_PC)
CI_res <- rbind(ubr1_CI, ubr2_CI, gray_CI, gcsi_CI, gdsc_CI, ctrp_CI, ccle_CI)

# save results
save(PC_res, CI_res,
     file = "DrugResponse/results/data/ARCHE_associations.RData")

###########################################################
# Identify Class A Biomarkers
###########################################################

# Class A biomarkers: abs(PCC > 0.6) & FDR < 0.05 in 1 PSet
ClassA <- PC_res[which((abs(PC_res$pc >= 0.6)) & PC_res$FDR < 0.05),]
ClassA <- ClassA[order(ClassA$pc, decreasing = T),]
ClassA$rank <- factor(1:nrow(ClassA), levels = c(1:nrow(ClassA)))


###########################################################
# Check and remove discordant associations in other PSets
###########################################################

##### BASED ON ALL CLASS A BIOMARKERS BEING (+)ve PCC

# Class A associations across PSets
toPlot <- PC_res[PC_res$pairs %in% ClassA$pairs,]
keep <- names(table(toPlot$pair)[table(toPlot$pair)>1])
toPlot <- toPlot[toPlot$pair %in% keep,]

# plot discordant associations
plot_ClassA_allAssociations(toPlot)

# remove biomarkers with discordant associations
ClassA <- ClassA[-which(ClassA$pair %in% toPlot$pairs[toPlot$pc < 0]),]
toPlot <- toPlot[toPlot$pairs %in% ClassA$pairs,]

# save Class A biomarkers
write.csv(ClassA, file = "DrugResponse/results/data/ClassA_Biomarkers.csv", quote = F, row.names = F)
write.csv(toPlot, file = "DrugResponse/results/data/ClassA_allAssociations.csv", quote = F, row.names = F)

###########################################################
# Plots for Class A biomarkers
###########################################################

# plot Class A biomarker associations
plot_ClassA_biomarkersAssociations(ClassA)

# plot Class A biomarker associations across PSets
plot_ClassA_associationsAcrossPSets(toPlot)

# plot indiv scatter plots for Class A biomarker associations across PSets
for (pair in ClassA$pair) {
    plot_ClassA_indivPlot(pair)
}
