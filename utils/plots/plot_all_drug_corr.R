#' Creates dataframe for drug correlation and plotting
#'
#' @param pset1 dataframe. Drug sensitivity dataframe1
#' @param pset2 dataframe. Drug sensitivity dataframe2
#' @return dataframe object for correlation and plotting.
#' 
format_drug_pset <- function(pset1, pset2) {

    # get common samples
    common_cell <- intersect(colnames(pset1), colnames(pset2))
    common_drug <- intersect(rownames(pset1), rownames(pset2))
    pset1 <- pset1[match(common_drug, rownames(pset1)),match(common_cell, colnames(pset1))]
    pset2 <- pset2[match(common_drug, rownames(pset2)),match(common_cell, colnames(pset2))]

    # set labels and melt
    pset1$drug <- rownames(pset1)
    pset2$drug <- rownames(pset2)
    pset1 <- melt(pset1)
    pset2 <- melt(pset2)
    pset1$pairs <- paste0(pset1$variable, "_", pset1$drug)
    pset2$pairs <- paste0(pset2$variable, "_", pset2$drug)

    # # scatter plot of drug response difference for CTRP and GDSC 
    toPlot <- data.frame(
        pair = pset1$pairs, 
        pset1 = pset1$value, 
        pset2 = pset2$value[match(pset1$pairs, pset2$pairs)]
    )
    
    return(toPlot)
}

#' Compile and plot all drug response correlation plots for PSet pairs
#' 
#' @param ____sen dataframe. Sensitivity data for each PSet
#' @param label string. "All" for all drugs, "BCa" for BCa relevant drugs only
#' @param corr_res dataframe. To store correlations
#' 
plot_all_drug_corr <- function(ubr1_sen, ubr2_sen, gray_sen, gcsi_sen, gdsc_sen, ctrp_sen, ccle_sen, label, corr_res) {

    # all correlation plots
    p1 <- plot_drug_corr(ubr1_sen, ubr2_sen, "UBR1", "UBR2", corr_res) 
    p2 <- plot_drug_corr(ubr1_sen, gray_sen, "UBR1", "GRAY", corr_res) 
    p3 <- plot_drug_corr(ubr1_sen, gcsi_sen, "UBR1", "gCSI", corr_res) 
    p4 <- plot_drug_corr(ubr1_sen, gdsc_sen, "UBR1", "GDSC2", corr_res) 
    p5 <- plot_drug_corr(ubr1_sen, ctrp_sen, "UBR1", "CTRP", corr_res) 
    p6 <- plot_drug_corr(ubr1_sen, ccle_sen, "UBR1", "CCLE", corr_res) 
    p7 <- plot_drug_corr(ubr2_sen, gray_sen, "UBR2", "GRAY", corr_res) 
    p8 <- plot_drug_corr(ubr2_sen, gcsi_sen, "UBR2", "gCSI", corr_res) 
    p9 <- plot_drug_corr(ubr2_sen, gdsc_sen, "UBR2", "GDSC2", corr_res) 
    p10 <- plot_drug_corr(ubr2_sen, ctrp_sen, "UBR2", "CTRP", corr_res) 
    p11 <- plot_drug_corr(ubr2_sen, ccle_sen, "UBR2", "CCLE", corr_res)
    p12 <- plot_drug_corr(gray_sen, gcsi_sen, "GRAY", "gCSI", corr_res) 
    p13 <- plot_drug_corr(gray_sen, gdsc_sen, "GRAY", "GDSC2", corr_res) 
    p14 <- plot_drug_corr(gray_sen, ctrp_sen, "GRAY", "CTRP", corr_res) 
    p15 <- plot_drug_corr(gray_sen, ccle_sen, "GRAY", "CCLE", corr_res) 
    p16 <- plot_drug_corr(gcsi_sen, gdsc_sen, "gCSI", "GDSC2", corr_res) 
    p17 <- plot_drug_corr(gcsi_sen, ctrp_sen, "gCSI", "CTRP", corr_res) 
    p18 <- plot_drug_corr(gcsi_sen, ccle_sen, "gCSI", "CCLE", corr_res) 
    p19 <- plot_drug_corr(gdsc_sen, ctrp_sen, "GDSC2", "CTRP", corr_res) 
    p20 <- plot_drug_corr(gdsc_sen, ccle_sen, "GDSC2", "CCLE", corr_res) 
    p21 <- plot_drug_corr(ctrp_sen, ccle_sen, "CTRP", "CCLE", corr_res)

    # output plots
    filename <- paste0("data/results/figures/3-DataExploration/drug_corr_", label, "Drugs.png")
    png(filename, width = 15, height = 13, res = 600, units = "in")
    print(
        grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21,
                ncol = 6, nrow = 6,
                layout_matrix = rbind(c(1, NA, NA, NA, NA, NA),
                                    c(2,  7, NA, NA, NA, NA),
                                    c(3,  8, 12, NA, NA, NA),
                                    c(4,  9, 13, 16, NA, NA),
                                    c(5, 10, 14, 17, 19, NA),
                                    c(6, 11, 15, 18, 20, 21)))
    )
    dev.off()

    # return the correlations
    return(corr_res)
}