#' Plot distribution of ARCHE scores
#' 
plot_ARCHE_scores <- function(toPlot) {

    png("DataExploration/results/figures/ARCHE_scores.png", width = 9, height = 4, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = Value)) + geom_density(fill = random_lightblue) + 
            facet_grid(.~ARCHE) + xlim(c(-255, 255)) +
            theme_classic() + geom_vline(xintercept = 0, linetype = 'dashed', color = "gray") +
            theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
            labs(y = 'Density', x = 'ARCHE Score'))
    dev.off()
}

#' Plot overlapping BCa CCLs in psets
#' 
plot_BCa_CCLs_overlap <- function(toPlot) {

    png("DataExploration/results/figures/BCa_CCLs_overlap.png", width = 4, height = 8, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = variable, y = reorder(sample, -as.numeric(factor(sample))), fill = factor(value))) + 
        geom_tile(color = "black") + 
        scale_fill_manual(values = c("0" = "white", "1" = "#ABC8C0"), labels = c("False", "True")) +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = "PSet", y = "Breast Cancer Cell Line", fill = "In PSet"))
    dev.off()
}

#' Plot overlapping BCa drugs in psets
#' 
plot_BCa_drugs_overlap <- function(toPlot) {

    png("DataExploration/results/figures/BCa_drugs_overlap.png", width=6, height=6, units='in', res = 600, pointsize=80)
    print(ggplot(toPlot, aes(x = PSet, y = Drug, fill = Present)) +
    geom_tile(color = "black") +
    scale_fill_manual(values = c(random_blue, "white")) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
            panel.grid = element_blank(),
            legend.key.size = unit(0.5, 'cm')) +
    labs(fill = "Drug\nPresent"))
    dev.off()

}

#' Plot correlation of drug AAC for pset pairs
#' 
#' @param pset1 dataframe. Drug sensitivity dataframe1
#' @param pset2 dataframe. Drug sensitivity dataframe2
#' @param label1 string. PSet1 label
#' @param label2 string. PSet2 label
#' @param corr_res dataframe. To store correlations
#' @return ggplot object of correlation plot
#' 
plot_drug_corr <- function(pset1, pset2, label1, label2, corr_res) {
    
    message(paste("***Working on", label1, "and", label2))
    toPlot <- format_drug_pset(pset1, pset2)

    # correlation and store results
    corr <- cor(toPlot$pset1, toPlot$pset2,  method = "pearson", use = "complete.obs")
    corr_res <- rbind(corr_res,
                      data.frame(PSet1 = label1, PSet2 = label2, PCC = corr))

    # create corr plot
    p <- ggplot(toPlot, aes(x = pset1, y = pset2)) + 
            geom_smooth(method=lm, show.legend = FALSE, color = "#046C9A") + 
            geom_point(shape = 21, size = 2.5, color = "black", fill = "#899DA4") + 
            geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") + 
            xlim(c(0, 1)) + ylim(c(0, 1)) + 
            theme_classic() + 
            theme(legend.key.size = unit(0.5, 'cm')) +
            geom_text(x = 0.2, y = 1, label = paste("PCC:", round(corr, digits = 3)), color = "black") +
            labs(x = label1, y = label2)
    return(p)
}

#' Compile and plot all drug response correlation plots for PSet pairs
#' 
#' @param ____sen dataframe. Sensitivity data for each PSet
#' @param label string. "All" for all drugs, "BCa" for BCa relevant drugs only
#' @param corr_res dataframe. To store correlations
#' 
plot_allDrugCorr <- function(ubr1_sen, ubr2_sen, gray_sen, gcsi_sen, gdsc_sen, ctrp_sen, ccle_sen, label, corr_res) {

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
    filename <- paste0("DataExploration/results/figures/drug_corr_", label, "Drugs.png")
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

#' Plot correlation results from RNA-Seq pset analysis
#' 
plot_rna_corr <- function(p1, p2, p3, p4, p5, p6) {
    # initialize correlation matrix
    psets <- c("UBR2", "GRAY", "gCSI", "CCLE")
    mat <- matrix(NA, nrow = 4, ncol = 4, dimnames = list(psets, psets))

    # fill in matrix values
    mat[1,2] <- mat[2,1] <- p1
    mat[1,3] <- mat[3,1] <- p2
    mat[1,4] <- mat[4,1] <- p3
    mat[2,3] <- mat[3,2] <- p4
    mat[2,4] <- mat[4,2] <- p5
    mat[3,4] <- mat[4,3] <- p6

    # set diagonal as 1
    diag(mat) <- 1

    # set palette for plotting
    col_fun <- colorRamp2(c(-1, 0, 1), c(binary_pal[2], "white", binary_pal[1]))

    # function to create heatmap
    plot_heatmap <- function(mat, legend_label) {
        plot <- Heatmap(mat,
            col = col_fun,
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            heatmap_legend_param = list(
                title = legend_label,
                color_bar = "continuous"
            ),
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(round(mat[i, j], 2), x, y, gp = gpar(fontsize = 8))
            }
        )
        return(plot)
    }

    png("DataExploration/results/figures/rna_corr.png", width=125, height=100, units='mm', res = 600, pointsize=80)
    print(plot_heatmap(mat, "Spearman"))
    dev.off()
}