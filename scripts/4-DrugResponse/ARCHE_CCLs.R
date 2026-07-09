# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(matrixStats)
    library(meta)
    library(ggh4x)
    library(reshape2)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(tidyverse)
    library(readxl)
    library(data.table)
    library(patchwork)
    library(ggnewscale)
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/plots/drug_response_ccls.R")
source("utils/palettes.R")
source("utils/bca_drugs.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/ccl_benchmarks.R")
source("utils/plots/drug_response_pdx.R")

###########################################################
# Prepare metadata
###########################################################

# read in sample metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]

# get cell lines
c_meta <- meta[meta$type == "cell_line", ]

###########################################################
# Load in cell line data
###########################################################

# zscores
zscore_cells <- get_arche_scores(paste0("data/rawdata/all_scoring/cell_tcga.Zscore.txt"), c_meta)
normzs_cells <- znorm(zscore_cells)

# sumdevs
zscore_cells_sumdev <- get_arche_sumdevs(zscore_cells, "zscore_cells", plot = TRUE)
normzs_cells_sumdev <- znorm(zscore_cells_sumdev)

# get drug sensitivity data
load("data/procdata/CCLs/sensitivity_data.RData")

###########################################################
# Compute PC of ARCHE-drug associations
###########################################################

# helper function to compute arche associations across all psets
arche_pc <- function(scores, label) {
    ubr1_PC <- compute_pc(scores, ubr1_sen, "UBR1")
    ubr2_PC <- compute_pc(scores, ubr2_sen, "UBR2")
    gray_PC <- compute_pc(scores, gray_sen, "GRAY")
    gcsi_PC <- compute_pc(scores, gcsi_sen, "gCSI")
    gdsc_PC <- compute_pc(scores, gdsc_sen, "GDSC2")
    ctrp_PC <- compute_pc(scores, ctrp_sen, "CTRP")
    ccle_PC <- compute_pc(scores, ccle_sen, "CCLE")

    # compile results
    PC_res <- rbind(ubr1_PC, ubr2_PC, gray_PC, gcsi_PC, gdsc_PC, ctrp_PC, ccle_PC)
    PC_res$Label <- label

    # map drugs
    PC_res <- map_drugs(PC_res)
    return(PC_res)
}

# zscores
pc_zscore_cells <- arche_pc(zscore_cells, "zscore")
pc_normzs_cells <- arche_pc(normzs_cells, "normzscr")

# subsetted zscores
pc_zscore_cells_sumdev <- arche_pc(zscore_cells_sumdev, "zscore_sumdev")
pc_normzs_cells_sumdev <- arche_pc(normzs_cells_sumdev, "normzscr_sumdev")

# save results
save(pc_zscore_cells, pc_normzs_cells, pc_zscore_cells_sumdev, pc_normzs_cells_sumdev,
     file = "data/results/data/4-DrugResponse/CCLs/ARCHE_CCLs_associations.RData")

###########################################################
# Identify Class A Biomarkers
###########################################################

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
    filepath <- paste0("data/results/data/4-DrugResponse/CCLs/ClassA_Biomarkers_", label, ".csv")
    write.csv(ClassA, file = filepath, quote = FALSE, row.names = FALSE)
    filepath <- paste0("data/results/data/4-Drugresponse/CCLs/ClassA_allAssociations_", label, ".csv")
    write.csv(toPlot, file = filepath, quote = FALSE, row.names = FALSE)

    # plot Class A heatmap (drug in >1 PSet)
    plot_ClassA_heatmap(toPlot, "Multi", label)

    # plot Class A associations in 1 PSet
    toPlot <- ClassA[-which(ClassA$pairs %in% toPlot$pairs),]
    toPlot <- map_drugs(toPlot)         # deal w drug names
    plot_ClassA_heatmap(toPlot, "Single", label)

}

get_classA(pc_zscore_cells, "zscore")
get_classA(pc_normzs_cells, "normzscr")

get_classA(pc_zscore_cells_sumdev, "zscore_sumdev")
get_classA(pc_normzs_cells_sumdev, "normzscr_sumdev")

###########################################################
# Identify Class B Biomarkers
###########################################################

# helper function to compile PCC and meta-estimates for ClassB biomarkers
get_classB <- function(PC_res, label) {

    # perform meta analysis and save results
    estimates <- compute_meta(PC_res)
    filepath <- paste0("data/results/data/4-DrugResponse/CCLs/meta_estimates", label, ".csv")
    write.csv(estimates, file = filepath, quote = FALSE, row.names = FALSE)

    #ClassB biomarkers: abs(TE > 0.25) & FDR < 0.05
    ClassB <- estimates[which(abs(estimates$TE) > 0.25 & estimates$FDR < 0.05),]
    ClassB <- ClassB[order(ClassB$TE, decreasing = TRUE),]
    ClassB$rank <- factor(1:nrow(ClassB), levels = c(1:nrow(ClassB)))
    
    keep <- PC_res[which(PC_res$pairs %in% ClassB$pair),]
    
    # combine PC_res and meta results
    toPlot <- data.frame(
        signature = c(keep$signature, ClassB$signature),
        drug = c(keep$drug, ClassB$drug),
        pairs = c(keep$pairs, ClassB$pair),
        estimate = c(keep$pc, ClassB$TE),
        upper = c(keep$upper, ClassB$upper),
        lower = c(keep$lower, ClassB$lower),
        FDR = c(keep$FDR, ClassB$FDR),
        pset = c(keep$pset, rep("Meta Estimate", nrow(ClassB)))
    )
    toPlot$meta <- factor(ifelse(toPlot$pset == "Meta Estimate", TRUE, FALSE), levels = c(TRUE, FALSE))
    toPlot$pset <- factor(toPlot$pset, levels = c(unique(PC_res$pset), "Meta Estimate"))

    # plot Class B biomarker associations as heatmap
    plot_ClassB_heatmap(toPlot, label)

    return(unique(toPlot$pairs))

}

zscore_cells_pairs <- get_classB(pc_zscore_cells, "zscore")
normzs_cells_pairs <- get_classB(pc_normzs_cells, "normzscr")

zscore_cells_sumdev_pairs <- get_classB(pc_zscore_cells_sumdev, "zscore_sumdev")
normzs_cells_sumdev_pairs <- get_classB(pc_normzs_cells_sumdev, "normzscr_sumdev")


###########################################################
# ARCHE2 vs ARCHE5
###########################################################

a2va5 <- compile_diff(zscore_cells_sumdev, "ARCHE2", "ARCHE5", a2va5_drug_pal, "A2vsA5")