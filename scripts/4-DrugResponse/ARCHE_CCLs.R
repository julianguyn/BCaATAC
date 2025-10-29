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

# load in BCa cell lines
samples <- get_cells()

# load in signature scores
signature_scores <- get_arche_cells()

# read in subtype information
true_subtype <- get_cell_subtype()

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
# Subset signature associations
###########################################################

# helper function to get PSet-specific signature scores
get_scores <- function(sen) {
    to_keep <- samples[which(samples$sample %in% colnames(sen)),]
    sig <- signature_scores[,which(colnames(signature_scores) %in% to_keep$sample)]
    sig <- sig[,order(colnames(sig))]
    return(sig)
}

# keep only CCLs with drug response per PSet
ubr1_sig <- get_scores(ubr1_sen)
ubr2_sig <- get_scores(ubr2_sen)
gray_sig <- get_scores(gray_sen)
gcsi_sig <- get_scores(gcsi_sen)
gdsc_sig <- get_scores(gdsc_sen)
ctrp_sig <- get_scores(ctrp_sen)
ccle_sig <- get_scores(ccle_sen)
tdxd_sig <- get_scores(tdxd_sen)

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

###########################################################
# Compute CI and PC of ARCHE-drug associations
###########################################################

# PC
ubr1_PC <- compute_pc(ubr1_sig, ubr1_sen, "UBR1")
ubr2_PC <- compute_pc(ubr2_sig, ubr2_sen, "UBR2")
gray_PC <- compute_pc(gray_sig, gray_sen, "GRAY")
gcsi_PC <- compute_pc(gcsi_sig, gcsi_sen, "gCSI")
gdsc_PC <- compute_pc(gdsc_sig, gdsc_sen, "GDSC2")
ctrp_PC <- compute_pc(ctrp_sig, ctrp_sen, "CTRP")
ccle_PC <- compute_pc(ccle_sig, ccle_sen, "CCLE")
tdxd_PC <- compute_pc(tdxd_sig, tdxd_sen, "TDXd")
erbb2_PC <- compute_pc(erbb2, erbb2_tdxd, "TDXd")

# CI (not used)
ubr1_CI <- compute_ci(ubr1_sig, ubr1_sen, "UBR1")
ubr2_CI <- compute_ci(ubr2_sig, ubr2_sen, "UBR2")
gray_CI <- compute_ci(gray_sig, gray_sen, "GRAY")
gcsi_CI <- compute_ci(gcsi_sig, gcsi_sen, "gCSI")
gdsc_CI <- compute_ci(gdsc_sig, gdsc_sen, "GDSC2")
ctrp_CI <- compute_ci(ctrp_sig, ctrp_sen, "CTRP")
ccle_CI <- compute_ci(ccle_sig, ccle_sen, "CCLE")

# compile results
PC_res <- rbind(ubr1_PC, ubr2_PC, gray_PC, gcsi_PC, gdsc_PC, ctrp_PC, ccle_PC)
CI_res <- rbind(ubr1_CI, ubr2_CI, gray_CI, gcsi_CI, gdsc_CI, ctrp_CI, ccle_CI)

# save results
save(PC_res, CI_res,
     file = "data/results/data/4-DrugResponse/ARCHE_CCLs_associations.RData")

###########################################################
# Identify Class A Biomarkers
###########################################################

# Class A biomarkers: abs(PCC > 0.5) & FDR < 0.05 in 1 PSet
ClassA <- PC_res[which((abs(PC_res$pc) >= 0.5) & PC_res$FDR < 0.05),]
ClassA <- ClassA[order(ClassA$pc, decreasing = T),]
ClassA$rank <- factor(1:nrow(ClassA), levels = c(1:nrow(ClassA)))


###########################################################
# Check and remove discordant associations in other PSets
###########################################################

# Class A associations across PSets
toPlot <- PC_res[PC_res$pairs %in% ClassA$pairs,]
keep <- names(table(toPlot$pair)[table(toPlot$pair)>1])
toPlot <- toPlot[toPlot$pair %in% keep,]

# plot all associations per ARCHE
plot_ClassA_allAssociations(toPlot, "ARCHE1", 8)
plot_ClassA_allAssociations(toPlot, "ARCHE2", 18)
plot_ClassA_allAssociations(toPlot, "ARCHE3", 20)
plot_ClassA_allAssociations(toPlot, "ARCHE4", 17)
plot_ClassA_allAssociations(toPlot, "ARCHE5", 22)
plot_ClassA_allAssociations(toPlot, "ARCHE6", 17)

# TODO: remove biomarkers with discordant associations?

# save Class A biomarkers
write.csv(ClassA, file = "data/results/data/4-DrugResponse/ClassA_Biomarkers.csv", quote = F, row.names = F)
write.csv(toPlot, file = "data/results/data/4-Drugresponse/ClassA_allAssociations.csv", quote = F, row.names = F)

###########################################################
# Plots for Class A biomarkers
###########################################################

# plot Class A heatmap (drug in >1 PSet)
plot_ClassA_heatmap(toPlot, "Multi")

# plot Class A associations in 1 PSet
toPlot <- ClassA[-which(ClassA$pairs %in% toPlot$pairs),]
toPlot <- map_drugs(toPlot)         # deal w drug names
plot_ClassA_heatmap(toPlot, "Single")

# plot Class A biomarker associations
ClassA <- map_drugs(ClassA)
plot_ClassA_biomarkersAssociations(ClassA)

# plot indiv scatter plots for Class A biomarker associations across PSets
#for (pair in ClassA$pair) {
#    plot_indivPlot(pair, "ClassA")     # rm for now, just plot what's needed
#}

###########################################################
# Identify Class B Biomarkers
###########################################################

# perform meta analysis and save results
estimates <- compute_meta(PC_res)
write.csv(estimates, file = "data/results/data/4-DrugResponse/meta_estimates.csv", quote = F, row.names = F)

#ClassB biomarkers: abs(TE > 0.4) & FDR < 0.05
ClassB <- estimates[which(abs(estimates$TE) > 0.4 & estimates$FDR < 0.05),]
ClassB <- ClassB[order(ClassB$TE, decreasing = T),]
ClassB$rank <- factor(1:nrow(ClassB), levels = c(1:nrow(ClassB)))


###########################################################
# Plots for Class B biomarkers
###########################################################

# helper function to compile PCC and meta-estimates for ClassB biomarkers
compileClassB <- function(PC_res, ClassB) {
    
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

    return(toPlot)
}

# plot Class B biomarker associations
plot_ClassB_biomarkersAssociations(ClassB)

# plot Class B biomarker associations as forest plot
toPlot <- compileClassB(PC_res, ClassB)
plot_ClassB_forest(toPlot)

# plot Class B biomarker associations as heatmap
plot_ClassB_heatmap(toPlot)

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
