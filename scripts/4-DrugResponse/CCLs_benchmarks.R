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

###########################################################
# Load in data
###########################################################

# load ARCHE drug response associations
load("data/results/data/4-DrugResponse/CCLs/ARCHE_CCLs_associations.RData")

# load subtyping scores
load("data/results/data/3-DataExploration/ccls_subtyping_scores.RData")

# get drug sensitivity data
load("data/procdata/CCLs/sensitivity_data.RData")

###########################################################
# Format subtype scores
###########################################################

# helper function to format scores
format_scores <- function(scores) {
    scores$Subtype <- NULL
    scores <- t(scores) |> as.data.frame()
    return(scores)
}

ubr2_pam50 <- format_scores(ubr2_pam50)
ccle_pam50 <- format_scores(ccle_pam50)
gcsi_pam50 <- format_scores(gcsi_pam50)
gray_pam50 <- format_scores(gray_pam50)
ubr2_scmod2 <- format_scores(ubr2_scmod2)
ubr2_scmgene <- format_scores(ubr2_scmgene)

# check overlap with sensitivity dataframes
check_sample_overlap(colnames(ubr2_pam50), colnames(ubr2_sen), "pam50", "Sen")
check_sample_overlap(colnames(ccle_pam50), colnames(ctrp_sen), "pam50", "Sen")
check_sample_overlap(colnames(gray_pam50), colnames(gray_sen), "pam50", "Sen")
check_sample_overlap(colnames(gcsi_pam50), colnames(gcsi_sen), "pam50", "Sen")
check_sample_overlap(colnames(ubr2_scmod2), colnames(ubr2_sen), "scmod2", "Sen")
check_sample_overlap(colnames(ubr2_scmgene), colnames(ubr2_sen), "scmgene", "Sen")

###########################################################
# Compute PC of Subtype-drug associations
###########################################################

# compute subtype drug response associations across all psets
ubr2_pm_PC <- compute_pc(ubr2_pam50, ubr2_sen, "UBR2_pam50")
ubr2_sm_PC <- compute_pc(ubr2_scmod2, ubr2_sen, "UBR2_scmod2")
ubr2_sg_PC <- compute_pc(ubr2_scmgene, ubr2_sen, "UBR2_scmgene")
gcsi_pm_PC <- compute_pc(gcsi_pam50, gcsi_sen, "gCSI_pam50")
gray_pm_PC <- compute_pc(gray_pam50, gray_sen, "GRAY_pam50")
ccle_pm_PC <- compute_pc(ccle_pam50, ctrp_sen, "CCLE_pam50")

# compile results
PC_res <- rbind(ubr2_pm_PC, ubr2_sm_PC, ubr2_sg_PC, gcsi_pm_PC, gray_pm_PC, ccle_pm_PC)
save(PC_res, file = "data/results/data/4-DrugResponse/CCLs/subtypes_CCLs_associations.RData")


###########################################################
# Plot Basal / ER-/HER2- / ARCHE2&5 and paclitaxel
###########################################################

pairs <- c("Basal_Paclitaxel", "ER-/HER2-_Paclitaxel", "ARCHE2_Paclitaxel", "ARCHE5_Paclitaxel")

# helper function to subset PCC results
keep_pairs <- function(pc, label) {
    if (label == "subtype") {
        pc$label <- sub(".*_", "", pc$pset)
    } else {
        pc$label <- label
    }
    pc <- pc[pc$pairs %in% pairs,]
    return(pc)
}

# compile pc results
toPlot <- rbind(
    keep_pairs(PC_res, "subtype"),
    keep_pairs(pc_20k, "20k"),
    keep_pairs(pc_50k, "50k"),
    keep_pairs(pc_all, "all"),
    keep_pairs(pcnorm_20k, "20knorm"),
    keep_pairs(pcnorm_50k, "50knorm"),
    keep_pairs(pcnorm_all, "allnorm")
)
toPlot$sig <- ifelse(toPlot$FDR < 0.05, "FDR < 0.05", "FDR >= 0.05")
toPlot$pairs <- paste0(toPlot$pairs, "\n", toPlot$label)
toPlot$pset <- sub("_.*", "", toPlot$pset)

# plot PC associations
p <- ggplot(toPlot, aes(x = pset, y = pairs, fill = pc, size = -log(FDR), shape = sig)) +
    geom_point() +
    geom_text(data = toPlot, aes(label = sprintf("%.2f", pc)), color = "black", size = 2.5) +
    scale_shape_manual(values = c(21, 24)) +
    scale_size(range = c(2, 12)) +
    scale_fill_gradient2(
        low = "#BC4749",
        high = "#689CB0",
        mid = "#C2BBC9",
        limits = c(-1, 1)
    ) +
    theme_minimal() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 10)
    ) +
    labs(size = "-log(FDR)", y = "Subtype-Drug Pair", x = "", shape = "FDR\nSignificance", fill = "Pearson's\nCorrelation\nCoefficient")

filename <- paste0("data/results/figures/4-DrugResponse/benchmarks/paclitaxel.png")
png(filename, width=7, height=6, units='in', res = 600, pointsize=80)
p
dev.off()