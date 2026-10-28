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

    filename <- paste0("data/results/figures/2-MolecularSigAnalysis/deg/volcano_", label, "vsOther.png")
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

    filename <- paste0("data/results/figures/2-MolecularSigAnalysis/MYC/volcano_MYC_", label, "vsOther.png")
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

    filename <- paste0("data/results/figures/2-MolecularSigAnalysis/MYC/expression_per_",variable,"_", label, ".png")
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

#' Plot volcano plots with selected genes labeled
#' 
#' @param toPlot dataframe. Result from run_DEG()
#' @param genes string. Vector of genes to label
#' @param label string. Label for filename
#' @param subset boolean. If genes should be subsetted for DEG thresholds
#' 
plot_volcano_select_genes <- function(toPlot, genes, label, subset = FALSE) {

    toPlot <- toPlot[complete.cases(toPlot),]

    # get genes to label 
    toPlot$to_label <- ifelse(rownames(toPlot) %in% genes, "Yes", "No")

    # if subset == TRUE, only label abs(logFC) > 2 and padj < 5%
    if (subset == TRUE) {
        toPlot$to_label[toPlot$to_label == "Yes"] <- ifelse(toPlot$padj[toPlot$to_label == "Yes"] < 0.05, "Yes", "No")
        toPlot$to_label[toPlot$to_label == "Yes"] <- ifelse(abs(toPlot$log2FoldChange)[toPlot$to_label == "Yes"] > 1.5, "Yes", "No")
    }

    # label direction
    toPlot$direction <- ifelse(toPlot$log2FoldChange > 0, "Positive", "Negative")
    toPlot$direction <- factor(toPlot$direction, levels = c("Positive", "Negative"))

    filename <- paste0("data/results/figures/2-MolecularSigAnalysis/DMR/volcano_", label, "vsOther.png")
    png(filename, width = 6, height = 5, res = 600, units = "in")
    print(
        ggplot(toPlot, aes(x = log2FoldChange, y = -log(padj))) +
        geom_point(color = "grey") +
        geom_point(data = toPlot[toPlot$to_label == "Yes",], 
                    aes(x = log2FoldChange, y = -log(padj), color = direction)) +
        geom_text_repel(data = toPlot[toPlot$to_label == "Yes",],
                        aes(x = log2FoldChange, y = -log(padj), label = X),
                        size = 3.5, force = 10, max.overlaps = Inf, force_pull = 0.1) +
        geom_hline(yintercept = -log(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(1.5, -1.5), linetype = "dashed") +
        scale_color_manual(values = binary_pal) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste(label, "vs Other"))
    )
    dev.off()
}
#' Plot volcano plots highlighting downregulated genes
#' 
#' Volcano plot with significant genes coloured by direction, the number of
#' up/downregulated genes annotated, and selected genes labeled
#' @param deg dataframe. DEG results with gene names in column X
#' @param label string. Label for title and filename
#' @param genes string. Vector of genes to label
#' @param fc_thres int. Threshold for abs(logFC) significance
#' @param x_cap int. abs(logFC) beyond this is drawn at the axis edge
#' 
plot_volcano_down <- function(deg, label, genes, fc_thres = 2, x_cap = 10) {

    deg <- deg[complete.cases(deg$padj),]
    deg$dir <- "NS"
    deg$dir[deg$log2FoldChange > fc_thres & deg$padj < 0.05] <- "Up"
    deg$dir[deg$log2FoldChange < -fc_thres & deg$padj < 0.05] <- "Down"
    deg$dir <- factor(deg$dir, levels = c("Down", "Up", "NS"))

    # number of up/downregulated genes
    n_down <- sum(deg$dir == "Down")
    n_up <- sum(deg$dir == "Up")
    ymax <- max(-log10(deg$padj))

    to_label <- deg[deg$X %in% genes,]

    p <- ggplot(deg, aes(x = log2FoldChange, y = -log10(padj), color = dir)) +
        geom_point(data = deg[deg$dir == "NS",], size = 0.6) +
        geom_point(data = deg[deg$dir != "NS",], size = 0.9) +
        geom_point(data = to_label, size = 1.8, shape = 21, fill = NA, color = "black", stroke = 0.6) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
        geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed", linewidth = 0.3) +
        geom_label_repel(data = to_label, aes(label = X), color = "black", size = 3,
                        fill = alpha("white", 0.7), linewidth = 0, label.padding = unit(0.1, "lines"),
                        max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3, box.padding = 0.4) +
        annotate("text", x = -Inf, y = Inf, label = paste0(n_down, " down"), hjust = -0.2, vjust = 1.5,
                 size = 3.5, color = binary_pal3[2], fontface = "bold") +
        annotate("text", x = Inf, y = Inf, label = paste0(n_up, " up"), hjust = 1.2, vjust = 1.5,
                 size = 3.5, color = binary_pal3[1], fontface = "bold") +
        scale_color_manual(values = c("Down" = binary_pal3[2], "Up" = binary_pal3[1], "NS" = "grey85"), guide = "none") +
        scale_x_continuous(limits = c(-x_cap, x_cap), oob = scales::squish) + # extreme logFC from near-zero counts
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 11)) +
        labs(x = "log2 Fold Change", y = "-log10(padj)", title = paste(label, "vs Other"))

    filename <- paste0("data/results/figures/2-MolecularSigAnalysis/deg/volcano_down_", label, "vsOther.png")
    print(paste("Saved to", filename))
    ggsave(filename, p, width = 3.5, height = 3.5)
}
