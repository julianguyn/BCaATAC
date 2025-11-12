#' Plot maf summaries
#' 
#' Volcano plots with labeled top genes
#' @param arche string. ARCHE name
#' 
plot_mafSummary <- function(arche = "all") {

    if (arche == "all") {
        mafs <- merge_mafs(meta$SNV.File.Name)
    } else {
        mafs <- merge_mafs(meta$SNV.File.Name[meta$ARCHE == arche])
    }
    
    filename <- paste0("data/results/figures/2-MolecularSigAnalysis/maf/", arche, ".png")
    png(filename, width = 8, height = 6, res = 600, units = "in")
    print(
        plotmafSummary(maf = mafs, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    )
    dev.off()
}

#' Plot distribution of TMB by ARCHE
#' 
plot_tmb_boxplot <- function(tmb) {

    p <- ggplot(tmb, aes(x = ARCHE, y = total_perMB_log, fill = ARCHE)) +
    geom_boxplot() + geom_jitter(width = 0.1, alpha = 0.5) +
    geom_signif(
            y_position = 0.7,
            xmin = 3, xmax = 6,
            annotation = round(sig_tukey$'p adj', 4),
            tip_length = 0.01,
            textsize = 4) +
    theme_classic() + 
    theme(
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "none") +
    scale_fill_manual(NULL, values = ARCHE_pal) +
    labs(y = "TMB/MB (log10)", x = "Assigned ARCHE")

    png(paste0("data/results/figures/2-MolecularSigAnalysis/maf/tmb_boxplot.png"), width = 4, height = 4, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Waterfall plot ranked by TMB coloured by ARCHE
#' 
# plot rank TMB and ARCHE
plot_tmb_waterfall <- function(tmb) {
  tmb$Rank <- 1:nrow(tmb) |> as.factor()

  p <- ggplot(tmb, aes(x = Rank, y = total_perMB_log, fill = ARCHE)) +
    geom_col(color = "black", linewidth = 0.3) +
    scale_fill_manual(values = ARCHE_pal) +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(y = "TMB/MB (log10)", x = "Tumour Sample")

  png(paste0("data/results/figures/2-MolecularSigAnalysis/maf/tmb_waterfall.png"), width = 6, height = 4, res = 600, units = "in")
  print(p)
  dev.off()
}

#' Plot detection of BCa-relevant mutations
#' 
plot_BCa_mutations <- function(toPlot) {
    
    # save labels for plotting
    toPlot <- reshape2::melt(toPlot)
    toPlot$ARCHE <- meta$ARCHE[match(toPlot$Var2, meta$snv_label)]
    toPlot$Subtype <- meta$Subtype[match(toPlot$Var2, meta$snv_label)]
    toPlot$rank <- meta$rank[match(toPlot$Var2, meta$snv_label)]

    # format for plotting
    toPlot <- toPlot[order(toPlot$rank),]
    toPlot$rank <- factor(toPlot$rank, levels = c(75:1))
    toPlot$value <- ifelse(toPlot$value >= 1, 1, 0)
    toPlot$value <- factor(toPlot$value, levels = c(1, 0))

    # format annotation bars
    toPlot_map <- toPlot[,colnames(toPlot) %in% c("Var2", "ARCHE", "Subtype", "rank")]
    toPlot_map <- toPlot_map[order(toPlot_map$rank),]
    toPlot_map$rank <- factor(toPlot_map$rank, levels = c(75:1))

    # plot heatmap
    p1 <- ggplot(toPlot) + geom_tile(aes(x = Var1, y = rank, fill = value), color = "gray") +
    scale_fill_manual("Mutation Status", values = c(random_blue, "white"), labels = c("Mutated", "Not Mutated")) +
    geom_hline(yintercept = c(16.5, 24.5, 40.5, 49.5, 58.5), color = "black") +
    theme_void() + 
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1), 
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
    labs(x = "Genes", y = "Tumour Sample\n")

    # signature annotation bar
    p2 <- ggplot(toPlot_map) + geom_tile(aes(x = 1, y = rank, fill = ARCHE)) +
    theme_void() +
    scale_fill_manual(values = ARCHE_pal) +
    theme(axis.title.x = element_text(size=12, angle = 90, vjust = 0.5)) + 
    labs(fill = "ARCHE            ", x = "              ")

    # subtype annotation bar
    p3 <- ggplot(toPlot_map) + geom_tile(aes(x = 1, y = rank, fill = Subtype)) +
    theme_void() +
    scale_fill_manual(values = subtype_pal) +
    theme(axis.title.x = element_text(size=12, angle = 90, vjust = 0.5)) + 
    labs(fill = "Subtype          ", x = "              ")

    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    l3 <- as_ggplot(get_legend(p3))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")
    p3 <- p3+theme(legend.position = "none")

    png("data/results/figures/2-MolecularSigAnalysis/maf/BCa_mutations.png", width = 6, height = 6, res = 600, units = "in")
    print(
        grid.arrange(p1, p2, p3, l1, l2, l3, ncol = 19, nrow = 6,
                layout_matrix = rbind(c(1,1,1,1,1,1,1,1,1,1,1,2,3,4,4,4,4), 
                                    c(1,1,1,1,1,1,1,1,1,1,1,2,3,5,5,5,5), 
                                    c(1,1,1,1,1,1,1,1,1,1,1,2,3,5,5,5,5), 
                                    c(1,1,1,1,1,1,1,1,1,1,1,2,3,6,6,6,6),
                                    c(1,1,1,1,1,1,1,1,1,1,1,2,3,6,6,6,6),
                                    c(1,1,1,1,1,1,1,1,1,1,1,2,3,NA,NA,NA,NA)))
    )
    dev.off()
}
