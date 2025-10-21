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
