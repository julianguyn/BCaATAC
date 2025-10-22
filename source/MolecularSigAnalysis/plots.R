#' Plot maf summaries
#' 
#' Volcano plots with labeled top genes
#' @param arche string. ARCHE name
#' 
plot_mafSummary <- function(arche = "all") {

    if (arche == "all") {
        mafs <- merge_mafs(meta$SNV.File.Name)
    } else {
        mafs <- merge_mafs(meta$SNV.File.Name[meta$Signature == arche])
    }
    
    filename <- paste0("MolecularSigAnalysis/results/figures/maf/", arche, ".png")
    png(filename, width = 8, height = 6, res = 600, units = "in")
    print(
        plotmafSummary(maf = mafs, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    )
    dev.off()
}

#' Plot detection of BCa-relevant mutations
#' 
plot_BCa_mutations <- function(toPlot) {
    
    # save labels for plotting
    toPlot <- reshape2::melt(toPlot)
    toPlot$Signature <- meta$Signature[match(toPlot$Var2, meta$snv_label)]
    toPlot$Subtype <- meta$Subtype[match(toPlot$Var2, meta$snv_label)]
    toPlot$rank <- meta$rank[match(toPlot$Var2, meta$snv_label)]

    # format for plotting
    toPlot <- toPlot[order(toPlot$rank),]
    toPlot$rank <- factor(toPlot$rank, levels = c(75:1))
    toPlot$value <- ifelse(toPlot$value >= 1, 1, 0)
    toPlot$value <- factor(toPlot$value, levels = c(1, 0))

    # format annotation bars
    toPlot_map <- toPlot[,colnames(toPlot) %in% c("Var2", "Signature", "Subtype", "rank")]
    toPlot_map <- toPlot_map[order(toPlot_map$rank),]
    toPlot_map$rank <- factor(toPlot_map$rank, levels = c(75:1))

    # plot heatmap
    p1 <- ggplot(toPlot) + geom_tile(aes(x = Var1, y = rank, fill = value), color = "gray") +
    scale_fill_manual("Mutation Status", values = c("#533B4D", "white"), labels = c("Mutated", "Not Mutated")) +
    geom_hline(yintercept = c(16.5, 24.5, 40.5, 49.5, 58.5), color = "black") +
    theme_void() + 
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1), 
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
    labs(x = "Genes", y = "Tumour Sample\n")

    # signature annotation bar
    p2 <- ggplot(toPlot_map) + geom_tile(aes(x = 1, y = rank, fill = Signature)) +
    theme_void() +
    scale_fill_manual(values = ARCHE_pal) +
    theme(axis.title.x = element_text(size=12, angle = 90, vjust = 0.5)) + 
    labs(fill = "Signature          ", x = "              ")

    # subtype annotation bar
    p3 <- ggplot(toPlot_map) + geom_tile(aes(x = 1, y = rank, fill = Subtype)) +
    theme_void() +
    scale_fill_manual(values = subtype_pal) +
    theme(axis.title.x = element_text(size=12, angle = 90, vjust = 0.5)) + 
    labs(fill = "Subtype      ", x = "              ")

    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    l3 <- as_ggplot(get_legend(p3))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")
    p3 <- p3+theme(legend.position = "none")

    png("MolecularSigAnalysis/results/figures/maf/BCa_mutations.png", width = 6, height = 6, res = 600, units = "in")
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

#' Plot pheatmaps of mutational signatures
#' 
plot_pheatmap <- function(mat, label) {
    
    png(paste0("MolecularSigAnalysis/results/figures/molecularsig/pheatmap_", label, ".png"), width = 10, height = 10, res = 600, units = "in")
    print(
        pheatmap::pheatmap(mat = mat, 
                   cluster_rows = FALSE, 
                   main = paste("Cosine Similarity against", label, "Mutational Signatures"))
    )
    dev.off()
}

#' Plot heatmaps of mutational signatures
#' 
plot_heatmap_mutsig <- function(toPlot, type) {

    # add atac
    toPlot <- format_sig(toPlot)
  
    x_lab <- ifelse(type == "og30", "                     ", "                ")
    
    p1 <- ggplot(toPlot) + geom_tile(aes(x = variable, y = rank, fill = value), color = "gray") +
        scale_fill_distiller("Cosine Similarity", palette = "RdYlBu") +
        geom_hline(yintercept = c(16.5, 24.5, 40.5, 49.5, 58.5), linetype = "dotted", color = "black") +
        theme_void() + 
        theme(axis.text.x = element_text(size =8, angle = 90, vjust = 0.5, hjust=1), 
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
        labs(x = "\nMutational Signatures", y = "Tumour Sample\n")
    
    # signature annotation bar
    p2 <- ggplot(toPlot) + geom_tile(aes(x = 1, y = rank, fill = Signature)) +
        theme_void() +
        scale_fill_manual(values = ARCHE_pal) +
        theme(axis.title.x = element_text(size=12, angle = 90, vjust = 0.5)) + 
        labs(fill = "Signature          ", x = x_lab)
    
    # subtype annotation bar
    p3 <- ggplot(toPlot) + geom_tile(aes(x = 1, y = rank, fill = Subtype)) +
        theme_void() +
        scale_fill_manual(values = subtype_pal) +
        theme(axis.title.x = element_text(size=12, angle = 90, vjust = 0.5)) + 
        labs(fill = "Subtype      ", x = x_lab)
    
    
    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    l3 <- as_ggplot(get_legend(p3))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")
    p3 <- p3+theme(legend.position = "none")
    
    png(paste0("MolecularSigAnalysis/results/figures/molecularsig/heatmap_", type, ".png"), width = 8, height = 6, res = 600, units = "in")
    print(grid.arrange(p1, p2, p3, l1, l2, l3, ncol = 25, nrow = 6,
            layout_matrix = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,4,4,4), 
                                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,4,4,4), 
                                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,5,5,5,5), 
                                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,5,5,5,5),
                                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,6,6,6,6),
                                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,6,6,6,6))))
    dev.off()
}

#' Plot molecular sig - ARCHE correlations
#' 
plot_molecularsig_corr <- function(corr, label, type) {
  
    corr$ATAC.Sig <- factor(corr$ATAC.Sig, levels = c(paste0("ARCHE", 1:6)))
    
    p <- ggplot(corr) + geom_tile(aes(x = Mol.Sig, y = ATAC.Sig, fill = Corr), color = "gray") +
        scale_fill_gradient2("Spearman\nCorrelation",
                            low = "#AF4C5B", high = "#046C9A", mid = "white",
                            midpoint = 0, limits = c(-1, 1)) +
        theme_void() + 
        theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1), 
            axis.text.y = element_text(size = 8, hjust = 0.95),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
        labs(x = paste0(label, "\n"), y = "\nARCHE") + coord_flip()

    # set width and height
    if (type == "og30") {
        w <- 4
        h <- 7
    } else if (type == "v3") {
        w <- 3
        h <- 8.5
    } else {
        w <- 5
        h <- 8
    }
  
    png(paste0("MolecularSigAnalysis/results/figures/molecularsig/", type, "_corr.png"), width = w, height = h, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot boxplots of molsig correlations with ARCHEs
#' 
plot_corr_boxplots <- function(corr_v3, corr_hm) {

    # format toPlot
    corr_v3$Signature <- "Mutational Signatures"
    corr_hm$Signature <- "Hallmark Gene Set"
    toPlot <- rbind(corr_v3[,c(1,3,4)], corr_hm[,c(1,3,4)])
    toPlot$Signature <- factor(toPlot$Signature, levels = c("Mutational Signatures", "Hallmark Gene Set"))

    # get min and max values 
    mm <- toPlot %>%
    group_by(ATAC.Sig, Signature) %>%
    filter(Corr == max(Corr) | Corr == min(Corr)) %>%
    ungroup()

    png("MolecularSigAnalysis/results/figures/molecularsig/corr_boxplots.png", width = 6, height = 4.5, res = 600, units = "in")
    print(
        ggplot(toPlot, aes(x = ATAC.Sig, y = Corr, fill = Signature)) +
        geom_boxplot(color = "black") +
        theme_classic() + ylim(c(-1, 1)) + 
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
        scale_fill_manual(values = c("#077293", "#9AA5BD")) +
        geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed") +
        labs(y = "Spearman Correlation", x = "") +
        geom_point(data = mm, 
             aes(x = ATAC.Sig, y = Corr, group = Signature), 
             color = "#CE3333", size = 2,
             position = position_dodge(width = 0.75),
             show.legend = FALSE)
    )
    dev.off()
}

#' Plot volcano plots
#' 
#' Volcano plots with labeled top genes
#' @param deg dataframe. Result from run_DEG()
#' @param label string. Label for filename
#' @param fc_thres int. Threshold for abs(logFC) significance
#' @param n_top int. Number of top genes to label
#' 
plot_volcano <- function(deg, label, fc_thres = 2, n_top = 10) {

    deg$sig <- ifelse(abs(deg$log2FoldChange) > fc_thres & deg$padj < 0.05, TRUE, FALSE)

    # get n_top genes to label
    sig_genes <- deg[which(deg$sig == TRUE),]
    sig_genes <- sig_genes[order(abs(sig_genes$log2FoldChange), decreasing = TRUE),]
    top_genes <- sig_genes[1:n_top,]

    sig_genes$dir <- ifelse(sig_genes$log2FoldChange > 0, "Upregulated", "Downregulated")
    sig_genes$dir <- factor(sig_genes$dir, levels = c("Upregulated", "Downregulated"))

    filename <- paste0("MolecularSigAnalysis/results/figures/deg_volcanos/volcano_", label, "vsOther.png")
    png(filename, width = 6, height = 5, res = 600, units = "in")
    print(
        ggplot() + 
        geom_point(data = deg, aes(x = log2FoldChange, y = -log(padj)), color = "gray") +
        geom_point(data = sig_genes, aes(x = log2FoldChange, y = -log(padj), color = dir)) +
        geom_text_repel(data = top_genes, max.overlaps = 100, size = 2.5,
                        aes(x = log2FoldChange, y = -log(padj), label = rownames(top_genes))) +
        scale_color_manual("", values = binary_pal) +
        geom_hline(yintercept = -log(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
        #ylim(1:-0.15) + xlim(-1.55,1.55) + 
        theme_classic() + 
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste(label, "vs Other"))
    )
    dev.off()
}

#' Plot MYC volcano plots
#' 
#' Volcano plots with MYC labeled
#' @param toPlot dataframe. Result from run_DEG()
#' @param label string. Label for filename
#' 
plot_volcano_MYC <- function(toPlot, label) {

    toPlot$myc <- ifelse(rownames(toPlot) == "MYC", "MYC", "No")

    filename <- paste0("MolecularSigAnalysis/results/figures/MYC/volcano_MYC_", label, "vsOther.png")
    png(filename, width = 6, height = 5, res = 600, units = "in")
    print(
        ggplot(toPlot, aes(x = log2FoldChange, y = -log(padj))) +
        geom_point(color = "grey") +
        geom_point(data = toPlot[rownames(toPlot) == "MYC",], 
                    aes(x = log2FoldChange, y = -log(padj)), 
                    color = "red") +
        geom_text(data = toPlot[rownames(toPlot) == "MYC",],
                    aes(x = log2FoldChange, y = -log(padj), label = "MYC"),
                    vjust = -1, hjust = 0.5, size = 3.5) +
        geom_hline(yintercept = -log(0.05), linetype = "dashed") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste(label, "vs Other"))
    )
    dev.off()
}

#' MYC expression in tumours
#' 
#' Plot MYC expression by either ARCHE or subtype
#' @param toPlot dataframe.
#' @param variable string. "Signature" or "Subtype"
#' @param label string. Label for filename
#' 
plot_MYCexp <- function(toPlot, variable, label) {

    pals <- list(ARCHE_pal, subtype_pal)
    pal <- ifelse(variable == "Signature", pals[1], pals[2]) |> unlist()

    # set figure width
    w <- ifelse(label == "basal", 4, 6)

    filename <- paste0("MolecularSigAnalysis/results/figures/MYC/expression_per_",variable,"_", label, ".png")
    png(filename, width = w, height = 5, res = 600, units = "in")
    print(
        ggplot(myc, aes(x = .data[[variable]], y = MYC, fill = .data[[variable]])) +
            geom_boxplot() + 
            scale_fill_manual(values = pal) +
            theme_classic() +
            labs(y = "MYC expression")
    )
    dev.off()
}
