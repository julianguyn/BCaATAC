#' Plot distribution of ARCHE scores
#' 
plot_ARCHE_scores <- function(toPlot) {

    png("data/results/figures/3-DataExploration/ARCHE_scores.png", width = 9, height = 4, res = 600, units = "in")
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

    png("data/results/figures/3-DataExploration/BCa_CCLs_overlap.png", width = 4, height = 8, res = 600, units = "in")
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

    png("data/results/figures/3-DataExploration/BCa_drugs_overlap.png", width=6, height=6, units='in', res = 600, pointsize=80)
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

#' Plot correlation results from RNA-Seq pset analysis
#' 
plot_rna_corr <- function(p1, p2, p3, p4, p5, p6, corr) {
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

    filename <- paste0("data/results/figures/3-DataExploration/rna_", corr, ".png")
    png(filename, width=125, height=100, units='mm', res = 600, pointsize=80)
    print(plot_heatmap(mat, corr))
    dev.off()
}

#' Plot PAM50 subtypes
#' 
#' @param model string. Subtyping classification model for label
plot_bca_subtype <- function(toPlot, model) {

    # format dataframe for plotting
    toPlot$samples <- rownames(toPlot)
    toPlot <- reshape2::melt(toPlot, id.vars = "samples")

    # modify subtype pal for plotting
    subtype_pal <- c(
        "Basal" = "#AF4C5B", 
        "TNBC" = "#AF4C5B", 
        "ER-/HER2-" = "#AF4C5B", 
        "Her2" = "#EED4D3", 
        "LumA" = "#B3B4D0", 
        "ER+/HER2- Low Prolif" = "#B3B4D0",
        "LumB" = "#363E62", 
        "ER+/HER2- High Prolif" = "#363E62", 
        "Normal" = "#6365AF", 
        "Not Available" = "#eFeBF7")
    toPlot$value <- factor(toPlot$value, levels = names(subtype_pal))

    p1 <- ggplot(toPlot[toPlot$variable == "true_subtype",], aes(x = 1, y = samples, fill = value)) + 
        geom_tile(color = "white") + theme_void() +
        geom_text(aes(label = value), color = "black", size = 3) +
        scale_fill_manual(
            values = subtype_pal, 
            na.value = "#D1D7DD") +        
        theme(
                axis.text.x = element_blank(),
                axis.title.x = element_text(),
                axis.text.y = element_text(vjust = 0.5, hjust = 1),
                axis.title.y = element_text(angle = 90),
                strip.text.x = element_text(),
                legend.position = "none"
            ) +
        labs(x = "Subtype\n\n", y = "Cell Line", fill = "")
    
    p2 <- ggplot(toPlot[toPlot$variable != "true_subtype",], aes(x = variable, y = samples, fill = value)) + 
    geom_tile(color = "white") + theme_void() +
    geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_manual(values = subtype_pal, na.value = "#D1D7DD") +        
    theme(
            axis.text.x = element_text(vjust = 0.5),
            axis.title.x = element_text(),
            axis.text.y = element_blank(),
            axis.title.y = element_text(angle = 90),
            strip.text.x = element_text(),
            legend.key.size = unit(0.7, 'cm')
        ) +
    labs(x = "\nPSet", y = "", fill = "Subtype")

    filename <- paste0("data/results/figures/3-DataExploration/bcasubtypes_",model,"_true.png")
    png(filename, width = 8, height = 7, res = 600, units = "in")
    grid::grid.newpage()
    grid::grid.draw(
        grid.arrange(p1, p2, ncol = 4, nrow = 2, layout_matrix = rbind(c(1,2,2,2), c(1,2,2,2)))
        #ggarrange(p1, p2)
    )
    dev.off()
}
