# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(survcomp)
    library(wesanderson)
    library(ggplot2)
    library(ggh4x)
    library(reshape2)
    library(meta)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(dplyr)
    library(readxl)
    library(data.table)
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/plots/drug_response_ccls.R")
source("utils/palettes.R")
source("utils/bca_drugs.R")

###########################################################
# Load in data
###########################################################

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]

# remove komal dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "komal"), ]

c_meta <- meta[meta$type == "cell_line", ]
p_meta <- meta[meta$type == "PDX", ]

# load in arche scores
cells_20k <- get_arche_scores("cells", "k20", c_meta)
cells_50k <- get_arche_scores("cells", "k50", c_meta)
cells_all <- get_arche_scores("cells", "all", c_meta)

pdxs_20k <- get_arche_scores("pdxs", "k20", p_meta)
pdxs_50k <- get_arche_scores("pdxs", "k50", p_meta)
pdxs_all <- get_arche_scores("pdxs", "all", p_meta)

# get drug sensitivity data
load("data/procdata/CCLs/sensitivity_data.RData")

# get tdxd response data
tdxd <- get_tdxd()

###########################################################
# Format TDXd response data
###########################################################

# map cell lines names
rownames(tdxd) <- map_cells(tdxd$Cell.Line)
tdxd$Cell.Line <- NULL
tdxd <- t(tdxd) |> as.data.frame()

# keep common samples
tdxd <- tdxd[,which(colnames(tdxd) %in% samples$sample)]
tdxd_sen <- tdxd[,order(colnames(tdxd))]

###########################################################
# Load in gene count data (for tdxd analysis)
###########################################################

# get BCa gene counts
ubr2 <- get_pset_rna("UBR2")

# get gene counts metadata
gene_meta <- read.table("data/procdata/CCLs/rna/UBR2_RNA_meta.tsv", header = TRUE)

# keep only ERBB2
keep <- gene_meta$Ensembl[gene_meta$Gene.Symbol == "ERBB2"]
erbb2 <- ubr2[rownames(ubr2) == keep,,drop=FALSE] |> as.data.frame()

# keep common cells
common <- intersect(colnames(erbb2), colnames(tdxd_sen))
erbb2 <- erbb2[,match(common, colnames(erbb2))]
erbb2_tdxd <- tdxd_sen[,match(common, colnames(tdxd_sen))]

###########################################################
# Stratify HER2 and non-HER2 samples (for tdxd analysis)
###########################################################

# get samples
her2 <- true_subtype$sample[true_subtype$subtype == "Her2"]
nonher2 <- true_subtype$sample[true_subtype$subtype != "Her2"]

# subset signatures
her2 <- tdxd_sig[,colnames(tdxd_sig) %in% her2]
nonher2 <- tdxd_sig[,colnames(tdxd_sig) %in% nonher2]

# subset drug response
her2_tdxd <- tdxd_sen[,match(colnames(her2), colnames(tdxd_sen))]
nonher2_tdxd <- tdxd_sen[,match(colnames(nonher2), colnames(tdxd_sen))]

# todo:: implement tdxd associations with new arche scores
#tdxd_PC <- compute_pc(tdxd_sig, tdxd_sen, "TDXd")
#erbb2_PC <- compute_pc(erbb2, erbb2_tdxd, "TDXd")

###########################################################
# Compute PC of ARCHE-drug associations
###########################################################

# helper function to compute arche associations across all psets
arche_pc <- function(scores) {
    ubr1_PC <- compute_pc(scores, ubr1_sen, "UBR1")
    ubr2_PC <- compute_pc(scores, ubr2_sen, "UBR2")
    gray_PC <- compute_pc(scores, gray_sen, "GRAY")
    gcsi_PC <- compute_pc(scores, gcsi_sen, "gCSI")
    gdsc_PC <- compute_pc(scores, gdsc_sen, "GDSC2")
    ctrp_PC <- compute_pc(scores, ctrp_sen, "CTRP")
    ccle_PC <- compute_pc(scores, ccle_sen, "CCLE")

    # compile results
    PC_res <- rbind(ubr1_PC, ubr2_PC, gray_PC, gcsi_PC, gdsc_PC, ctrp_PC, ccle_PC)
    return(PC_res)
}

pc_20k <- arche_pc(cells_20k)
pc_50k <- arche_pc(cells_50k)
pc_all <- arche_pc(cells_all)

# save results
save(pc_20k, pc_50k, pc_all,
     file = "data/results/data/4-DrugResponse/ARCHE_CCLs_associations.RData")

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
    filepath <- paste0("data/results/data/4-DrugResponse/ClassA_Biomarkers_", label, ".csv")
    write.csv(ClassA, file = filepath, quote = FALSE, row.names = FALSE)
    filepath <- paste0("data/results/data/4-Drugresponse/ClassA_allAssociations_", label, ".csv")
    write.csv(toPlot, file = filepath, quote = FALSE, row.names = FALSE)

    # plot Class A heatmap (drug in >1 PSet)
    plot_ClassA_heatmap(toPlot, "Multi", label)

    # plot Class A associations in 1 PSet
    toPlot <- ClassA[-which(ClassA$pairs %in% toPlot$pairs),]
    toPlot <- map_drugs(toPlot)         # deal w drug names
    plot_ClassA_heatmap(toPlot, "Single", label)

    return(toPlot)
}

classA_20k <- get_classA(pc_20k, "20k")
classA_50k <- get_classA(pc_50k, "50k")
classA_all <- get_classA(pc_all, "all")

###########################################################
# Indiv plots for associations of interest
###########################################################

# ARCHE5 and Paclitaxel
plot_indivPlot("ARCHE5_Paclitaxel", cells_20k, "cells_20k")
plot_indivPlot("ARCHE5_Paclitaxel", cells_50k, "cells_50k")
plot_indivPlot("ARCHE5_Paclitaxel", cells_all, "cells_all")

# ARCHE1 and Olaparib
plot_indivPlot("ARCHE1_Olaparib", cells_20k, "cells_20k")
plot_indivPlot("ARCHE1_Olaparib", cells_50k, "cells_50k")
plot_indivPlot("ARCHE1_Olaparib", cells_all, "cells_all")

# ARCHE2 and Olaparib
plot_indivPlot("ARCHE2_Olaparib", cells_20k, "cells_20k")
plot_indivPlot("ARCHE2_Olaparib", cells_50k, "cells_50k")
plot_indivPlot("ARCHE2_Olaparib", cells_all, "cells_all")

# ARCHE2 and Topotecan
plot_indivPlot("ARCHE2_Topotecan", cells_20k, "cells_20k")
plot_indivPlot("ARCHE2_Topotecan", cells_50k, "cells_50k")
plot_indivPlot("ARCHE2_Topotecan", cells_all, "cells_all")

# ARCHE2 and SN-38
plot_indivPlot("ARCHE2_SN-38", cells_20k, "cells_20k")
plot_indivPlot("ARCHE2_SN-38", cells_50k, "cells_50k")
plot_indivPlot("ARCHE2_SN-38", cells_all, "cells_all")


###########################################################
# Identify Class B Biomarkers
###########################################################

# helper function to compile PCC and meta-estimates for ClassB biomarkers
get_classB <- function(PC_res, label) {

    # perform meta analysis and save results
    estimates <- compute_meta(PC_res)
    filepath <- paste0("data/results/data/4-DrugResponse/meta_estimates", label, ".csv")
    write.csv(estimates, file = filepath, quote = FALSE, row.names = FALSE)

    #ClassB biomarkers: abs(TE > 0.4) & FDR < 0.05
    ClassB <- estimates[which(abs(estimates$TE) > 0.4 & estimates$FDR < 0.05),]
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

    # plot Class B biomarker associations as forest plot
    #plot_ClassB_forest(toPlot, label)

    # plot Class B biomarker associations as heatmap
    plot_ClassB_heatmap(toPlot, label)

    return(toPlot)
}

classB_20k <- get_classB(pc_20k, "20k")
classB_50k <- get_classB(pc_50k, "50k")
classB_all <- get_classB(pc_all, "all")


###########################################################
# Combine TDXd data for plotting
###########################################################

# helper function to combine drug response and arche scores
combine_data <- function(signature_scores, tdxd) {
    combined <- t(rbind(tdxd, signature_scores)) |> as.data.frame()
    samples <- rownames(combined)
    combined <- sapply(combined, as.numeric) |> as.data.frame()
    rownames(combined) <- samples
    combined$subtype <- true_subtype$subtype[match(samples, true_subtype$sample)]
    return(combined)
}

all_c <- combine_data(tdxd_sig, tdxd)
her_c <- combine_data(her2, her2_tdxd)
non_c <- combine_data(nonher2, nonher2_tdxd)
erbb2_c <- combine_data(erbb2, erbb2_tdxd)

###########################################################
# Plot correlation of ARCHE and TDXd response
###########################################################

plot_tdxd_all(all_c, tdxd_PC, "Avg.IC50", "All Cells")
plot_tdxd_all(all_c, tdxd_PC, "Avg.IC50.Treps", "All Cells")
plot_tdxd_all(her_c, tdxd_PC, "Avg.IC50", "HER2 Cells")
plot_tdxd_all(her_c, tdxd_PC, "Avg.IC50.Treps", "HER2 Cells")
plot_tdxd_all(non_c, tdxd_PC, "Avg.IC50", "Non-HER2 Cells")
plot_tdxd_all(non_c, tdxd_PC, "Avg.IC50.Treps", "Non-Her2 Cells")

plot_tdxd_corr(erbb2_c, erbb2_PC, "ENSG00000141736", "Avg.IC50", save = TRUE)
plot_tdxd_corr(erbb2_c, erbb2_PC, "ENSG00000141736", "Avg.IC50.Treps", save = TRUE)
