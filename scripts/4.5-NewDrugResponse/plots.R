# test ARCHE scores with TCGA peaks as the bakground

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(matrixStats)
    library(ggh4x)
    library(reshape2)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(dplyr)
    library(readxl)
    library(data.table)
    library(patchwork)
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/palettes.R")
source("utils/bca_drugs.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/plots/drug_response_pdx.R")


# ---------------------------------------------------------
# Helper functions

# helper function to identify Class A biomarkers
get_classA <- function(PC_res, label) {
    # Class A biomarkers: abs(PCC > 0.5) & FDR < 0.05 in 1 PSet
    ClassA <- PC_res[which((abs(PC_res$pc) >= 0.5) & PC_res$FDR < 0.05),]
    ClassA <- ClassA[order(ClassA$pc, decreasing = T),]
    ClassA$rank <- factor(1:nrow(ClassA), levels = c(1:nrow(ClassA)))

    # get Class A associations across PSets
    toPlot <- PC_res[PC_res$pairs %in% ClassA$pairs,]
    keep <- names(table(toPlot$pair)[table(toPlot$pair)>1])
    toPlot <- toPlot[toPlot$pair %in% keep,]

    # save Class A biomarkers
    #filepath <- paste0("data/results/data/4-DrugResponse/CCLs/ClassA_Biomarkers_", label, ".csv")
    #write.csv(ClassA, file = filepath, quote = FALSE, row.names = FALSE)
    #filepath <- paste0("data/results/data/4-Drugresponse/CCLs/ClassA_allAssociations_", label, ".csv")
    #write.csv(toPlot, file = filepath, quote = FALSE, row.names = FALSE)

    # plot Class A heatmap (drug in >1 PSet)
    plot_ClassA_heatmap(toPlot, "test-Multi", label)

    # plot Class A associations in 1 PSet
    toPlot <- ClassA[-which(ClassA$pairs %in% toPlot$pairs),]
    toPlot <- map_drugs(toPlot)         # deal w drug names
    plot_ClassA_heatmap(toPlot, "test-Single", label)

    return(toPlot)
}


# helper function to plot individual plots
indiv_plots <- function(pair) {

    for (score in unique(cell_toPlot$Label)) {
        plot_indivPlot(pair, cell_toPlot[cell_toPlot$Label == score,], score, dir = "NewDrugResponse")
    }

}

###########################################################
# Load in cell line data
###########################################################

load("data/results/data/Misc/sample_tcga_ARCHE_drug_response_testing_all_scoring.RData")

# get drug sensitivity data
load("data/procdata/CCLs/sensitivity_data.RData")


###########################################################
# Indiv plots for associations of interest
###########################################################

# helper function to plot individual plots
indiv_plots <- function(pair) {

    for (score in unique(cell_toPlot$Label)) {
        plot_indivPlot(pair, cell_toPlot[cell_toPlot$Label == score,], score)
    }

}

# plot for pairs of interest
indiv_plots("ARCHE5_Paclitaxel")
