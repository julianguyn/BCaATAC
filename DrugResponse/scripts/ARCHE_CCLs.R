setwd("/Users/julianguyen/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(survcomp)
    library(wesanderson)
    library(ggplot2)
    library(ggh4x)
    library(reshape2)
    library(meta)
    library(ggh4x)
    library(ggpubr)
    library(grid)
    library(gridExtra)
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

# Class A biomarkers: abs(PCC > 0.65) & FDR < 0.05 in 1 PSet
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

# plot discordant associations
plot_ClassA_allAssociations(toPlot, "ARCHE1", 8)
plot_ClassA_allAssociations(toPlot, "ARCHE2", 18)
plot_ClassA_allAssociations(toPlot, "ARCHE3", 20)
plot_ClassA_allAssociations(toPlot, "ARCHE4", 17)
plot_ClassA_allAssociations(toPlot, "ARCHE5", 22)
plot_ClassA_allAssociations(toPlot, "ARCHE6", 17)

# TODO: remove biomarkers with discordant associations?

# save Class A biomarkers
write.csv(ClassA, file = "DrugResponse/results/data/ClassA_Biomarkers.csv", quote = F, row.names = F)
write.csv(toPlot, file = "DrugResponse/results/data/ClassA_allAssociations.csv", quote = F, row.names = F)

###########################################################
# Plots for Class A biomarkers
###########################################################

# make heatmap


plot_ClassA_heatmaps <- function(toPlot) {

    # subset signature and format plot
    toPlot <- toPlot[order(toPlot$pairs),]
    toPlot$pset <- factor(toPlot$pset, levels = names(PSet_pal))
    toPlot$pairs <- factor(toPlot$pairs, levels = rev(unique(toPlot$pairs)))

    # main heatmap (pset~pair)
    p1 <- ggplot(toPlot, aes(x = pset, y = pairs, fill = pc)) + 
        geom_tile(color = 'black') +
        geom_text(data = subset(toPlot, FDR < 0.05),
                aes(label = "*"), 
                vjust = 0.75, size = 4) +
        scale_fill_gradient2("Pearson's\nCorrelation\nCoefficient", 
                            low = binary_pal[2], 
                            high = binary_pal[1],
                            limits = c(-1, 1)) +
        scale_x_discrete(position = "top") +
        theme_void() +
        theme(
            axis.text.y = element_text(size=9, hjust=1, margin = margin(r = 3)), 
            axis.text.x = element_text(size = 9, margin = margin(b = 3)),
            legend.title = element_text(size = 9),
            axis.ticks = element_line(color = "gray", linewidth = 0.3),
            axis.ticks.length = unit(2, "pt")
            )

    # ARCHE annotation
    p2 <- ggplot(toPlot, aes(x = 1, y = pairs, fill = signature)) + 
        geom_tile(color = NA) +
        theme_void() + 
        scale_fill_manual("ARCHE", values = ARCHE_pal) +
        theme(axis.text.x = element_text(size = 9, margin = margin(b = 3))) + 
        labs(x = "")

    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")

    png("DrugResponse/results/figures/ClassA/heatmap.png", width = 8, height = 18, res = 600, units = "in")
    print(
        grid.arrange(p1, p2, l1, l2, ncol = 9, nrow = 2,
        layout_matrix = rbind(c(1,1,1,1,1,1,1,2,3), 
                              c(1,1,1,1,1,1,1,2,4)))
    )
    dev.off()
}

plot_ClassA_heatmaps(toPlot)


plot_ClassA_heatmaps <- function(toPlot, ARCHE, height) {

    # subset signature and format plot
    toPlot <- toPlot[toPlot$signature == ARCHE,]
    toPlot <- toPlot[order(toPlot$drug, decreasing = T),]
    toPlot$pset <- factor(toPlot$pset, levels = names(PSet_pal))
    toPlot$drug <- factor(toPlot$drug, levels = unique(toPlot$drug))
    filename <- paste0("DrugResponse/results/figures/ClassA/heatmap_", ARCHE, ".png")

    png(filename, width = 6, height = height, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = pset, y = drug, fill = pc)) + 
    geom_tile(color = 'black') +
    geom_text(data = subset(toPlot, FDR < 0.05),
            aes(label = "*"), 
            vjust = 0.75, size = 4) +
    scale_fill_gradient2("Pearson's\nCorrelation\nCoefficient", 
                         low = binary_pal[2], 
                         high = binary_pal[1],
                         limits = c(-1, 1)) +
    scale_x_discrete(position = "top") +
    theme_void() +
    theme(
        axis.text.y = element_text(size=9, hjust=1, margin = margin(r = 3)), 
        axis.text.x = element_text(size = 9, margin = margin(b = 3)),
        legend.title = element_text(size = 9),
        axis.ticks = element_line(color = "gray", linewidth = 0.3),
        axis.ticks.length = unit(2, "pt")))
    dev.off()
}

# plot heatmap
plot_ClassA_heatmaps(toPlot, "ARCHE1", 2)
plot_ClassA_heatmaps(toPlot, "ARCHE2", 5)
plot_ClassA_heatmaps(toPlot, "ARCHE3", 20)
plot_ClassA_heatmaps(toPlot, "ARCHE4", 17)
plot_ClassA_heatmaps(toPlot, "ARCHE5", 22)
plot_ClassA_heatmaps(toPlot, "ARCHE6", 17)

# plot Class A biomarker associations
plot_ClassA_biomarkersAssociations(ClassA)

# plot Class A biomarker associations across PSets
plot_ClassA_associationsAcrossPSets(toPlot)

# plot indiv scatter plots for Class A biomarker associations across PSets
for (pair in ClassA$pair) {
    plot_indivPlot(pair, "ClassA")
}

###########################################################
# Identify Class B Biomarkers
###########################################################

#ClassB biomarkers: abs(PCC > 0.4) & FDR < 0.05 in >1 PSet
to_keep <- PC_res[which((abs(PC_res$pc >= 0.4)) & PC_res$FDR < 0.05),]
to_keep <- names(table(to_keep$pairs)[table(to_keep$pairs)>1])

ClassB <- PC_res[which(PC_res$pairs %in% to_keep),]
ClassB <- ClassB[order(ClassB$pc, decreasing = T),]
ClassB$rank <- factor(1:nrow(ClassB), levels = c(1:nrow(ClassB)))


###########################################################
# Plots for Class B biomarkers
###########################################################

# plot Class B biomarker associations
plot_ClassB_biomarkersAssociations(ClassB)

# plot Class B biomarker associations across PSets
plot_ClassB_associationsAcrossPSets(ClassB)

# plot indiv scatter plots for Class A biomarker associations across PSets
for (pair in ClassB$pair) {
    plot_indivPlot(pair, "ClassB")
}

###########################################################
# Identify Class C Biomarkers
###########################################################

# perform meta analysis and save results
estimates <- run_meta(PC_res)
write.csv(estimates, file = "DrugResponse/results/data/meta_estimates.csv", quote = F, row.names = F)

#ClassC biomarkers: abs(TE > 0.4) & FDR < 0.05
ClassC <- estimates[which(abs(estimates$TE) > 0.4 & estimates$FDR < 0.05),]
ClassC <- ClassC[order(ClassC$TE, decreasing = T),]
ClassC$rank <- factor(1:nrow(ClassC), levels = c(1:nrow(ClassC)))


###########################################################
# Plots for Class C biomarkers
###########################################################

# plot Class C biomarker associations
plot_ClassC_biomarkersAssociations(ClassC)

# plot Class C biomarker associations as forest plot
plot_ClassC_forest(PC_res, ClassC)
