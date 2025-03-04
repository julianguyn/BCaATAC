# plot results from Sasha (top 10k sites in each CALS)

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(data.table)
})

###########################################################
# Load in data
###########################################################

cl1 <- fread("cfDNA/results/data/10k_sites/CAMA1_30/CAMA1_30.GC_corrected.coverage.tsv")
cl2 <- fread("cfDNA/results/data/10k_sites/MCF7_merged/MCF7_merged.GC_corrected.coverage.tsv")
pdo <- fread("cfDNA/results/data/10k_sites/BPTO95/BPTO95.GC_corrected.coverage.tsv")
xen <- fread("cfDNA/results/data/10k_sites/CAMA1_mouse_ctDNA_merged/CAMA1_mouse_ctDNA_merged.GC_corrected.coverage.tsv")


###########################################################
# Get CALS score
###########################################################

# function to get CALS score from Griffin output matrix
score_CALS <- function(griffin, label) {
    griffin <- data.frame(site = paste0("CALS-", 1:6),
                          sample = rep(label, 6),
                          #midpoint = griffin$'0',
                          #central = griffin$central_coverage,
                          CALS_score = 1-griffin$central_coverage)
    return(griffin)
}

cl1 <- score_CALS(cl1, "CAMA-1")
cl2 <- score_CALS(cl2, "MCF-7")
pdo <- score_CALS(pdo, "Organoid")
xen <- score_CALS(xen, "Xenograft")


###########################################################
# Merge tables for plotting
###########################################################

# average cell line coverage
ccl <- data.frame(site = paste0("CALS-", 1:6),
                  sample = rep("Cell Line", 6),
                  CALS_score = rowMeans(cbind(cl1$CALS_score, cl2$CALS_score)))

# merge results
toPlot <- rbind(ccl, pdo, xen)


###########################################################
# Plot coverage scores
###########################################################

# set palette
pal <- c("#ED7572", "#F4C5E1", "#A24848")

# bar plot of CALS score
png("cfDNA/results/figures/CALS_score_ER.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = site, y = CALS_score, fill = sample)) + 
    geom_bar(stat = "identity", position = position_dodge(), color = "black") + 
    theme_classic() + geom_hline(yintercept = 0) + 
    scale_fill_manual("Sample Type", values = pal) +
    theme(legend.key.size = unit(0.7, 'cm')) + 
    labs(x = "CALS", y = "CALS Score")
dev.off()