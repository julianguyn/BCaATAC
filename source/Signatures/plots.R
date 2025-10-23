#' ARCHE peak weight thresholds
#' 
#' Plot changes in peak weights from NMF for top n windows
#' @param peaks dataframe. Output of top_peaks()
#' @param arche string. ARCHE name
#' @param num_windows int. Number of windows kept
#' 
plot_peakWeights <- function(peaks, arche, num_windows) {

    n <- paste0(as.character(num_windows/1000), "k")

    # get start positions of plateaus
    diff <- find_delta0(peaks)

    filename <- paste0("Signatures/results/figures/peakWeights/",arche,"_", n, ".png")
    png(filename, width = 5, height = 5, res = 600, units = "in")
    print(
        ggplot(peaks, aes(x = rank, y = .data[[arche]])) +
        geom_line() +
        geom_point() +
        geom_vline(xintercept = diff$rank, linetype = "dashed") +
        geom_vline(xintercept = 10000) +
        theme_classic() +
        labs(x = "Peak Window Ranked by Weight", y = paste("Peak Window Weight for", arche))
    )
    dev.off()
}

#' ARCHE heatmap
#' 
#' Plot heatmap of ARCHE scores in TCGA BCa tumours
#' 
plot_ARCHE_heatmap <- function(mat) {

    p1 <- ggplot(mat, aes(x = variable, y = 1, fill = subtype)) + 
        geom_tile(color = NA) +
        theme_void() + 
        scale_fill_manual("BCa\nSubtype", values = subtype_pal,
                        labels = c("Basal", "Her2", "LuminalA", "LuminalB", "Normal", "NA")) +
        theme(axis.title.y = element_text(size=12)) + 
        labs(y = "              ")

    p2 <- ggplot(mat, aes(x = variable, y = 1, fill = signature_assign)) + 
        geom_tile(color = NA) +
        theme_void() + 
        scale_fill_manual("Assigned\nARCHE", values = ARCHE_pal) +
        theme(axis.title.y = element_text(size=12)) + 
        labs(y = "              ")

    p3 <- ggplot(mat, aes(x = variable, y = Signature, fill = value)) + 
        geom_tile(color = NA) +
        scale_fill_gradientn("ARCHE        \nExpression\nScore", colours = brewer.pal(9, "Blues")) + 
        theme_void() +
        theme(
            axis.text.y = element_text(size=11, margin = margin(r = 2)), 
            axis.title.x = element_text(size=12)
        ) + 
        labs(x = "Tumour Sample")

    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    l3 <- as_ggplot(get_legend(p3))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")
    p3 <- p3+theme(legend.position = "none")

    png("Signatures/results/figures/ARCHE_heatmap.png", width = 11, height = 4, res = 600, units = "in")
    print(
        grid.arrange(p1, p2, p3, l1, l2, l3, ncol = 9, nrow = 12,
        layout_matrix = rbind(c(1,1,1,1,1,1,1,4,5), c(2,2,2,2,2,2,2,4,5),
                            c(3,3,3,3,3,3,3,4,5), c(3,3,3,3,3,3,3,4,5),
                            c(3,3,3,3,3,3,3,4,5), c(3,3,3,3,3,3,3,4,5),
                            c(3,3,3,3,3,3,3,6,NA),c(3,3,3,3,3,3,3,6,NA),
                            c(3,3,3,3,3,3,3,6,NA),c(3,3,3,3,3,3,3,6,NA),
                            c(3,3,3,3,3,3,3,6,NA),c(3,3,3,3,3,3,3,NA,NA)))
    )
    dev.off()
}

#' Plot ARCHE window number
#' 
plot_ARCHE_peakInfo <- function(df) {

    df$ARCHE <- factor(df$ARCHE, levels = paste0("ARCHE", 6:1))
    filename <- "Signatures/results/figures/ARCHE_peakInfo.png"
    png(filename, width = 5, height = 5, res = 600, units = "in")
    print(
        ggplot(df, aes(x = num_windows, y = ARCHE)) +
            geom_bar(stat = "identity", fill = random_lightblue) +
            geom_text(data = df, aes(label = paste0(num_windows, " (", sum_peaks, "bp)"), x = 1.25e05)) +
            theme_classic() +
            labs(x = "Number of 500bp Windows")
    )
    dev.off()
}

#' ARCHE Upset plot
#' 
#' Plot Upset plot of overlapping peaks in ARCHEs
#' 
plot_ATAC_Upset <- function(m, label) {

    # set filter
    filter <- ifelse(label == "all", 2000000, 100000)
    m = m[comb_size(m) > filter]

    filename <- paste0("Signatures/results/figures/ATAC_Upset_", label, ".png")
    png(filename, width = 11, height = 5, res = 600, units = "in")
    print(
        UpSet(m, set_order = c(paste0("ARCHE", 1:6)), comb_order = order(-comb_size(m)),
        bg_col = "gray", comb_col = random_blue, 
        right_annotation = upset_right_annotation(m, gp = gpar(fill = random_blue)))
    )
    dev.off()
}

#' Plot annotatePeak results
#' 
plot_annotatePeak <- function(toPlot, label) {

    filename <- paste0("Signatures/results/figures/annotatePeak_", label, ".png")
    png(filename, width = 7, height = 5, res = 600, units = "in")
    print(
        ggplot(toPlot, aes(fill=Feature, y=Frequency, x=ARCHE)) + 
            geom_bar(position="fill", stat="identity", color = "black", size = 0.3) +
            scale_fill_manual(values = genfeat_pal) +
            theme_minimal() + 
            labs(y = "Percentage (%)")
        )
    dev.off()
}

#' Plot GREAT results
#' 
#' Plot lollipop graphs of GREAT enrcichment results
#' @param great dataframe. Output of runGREAT()
#' @param n_genes int. Number of annotated genes to filter by
#' @param label string. ARCHE label for filename
#' 
plot_GREAT <- function(great, n_genes, label) {

    # keep only top 30 with > n_genes
    great <- great[great$Total_Genes_Annotated > n_genes,]
    great <- great[order(great$Hyper_Fold_Enrichment, decreasing = TRUE),]
    toPlot <- great[1:30,]
    toPlot$name <- factor(toPlot$name, levels=rev(toPlot$name))

    p <- ggplot(toPlot, aes(x = name, y = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH)) +
        geom_point(aes(size = Total_Genes_Annotated), shape = 19) + 
        geom_segment(aes(x = name, xend = name, y = 0, yend = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH), size = 1) +
        coord_cartesian(clip = "off") + 
        coord_flip() +
        scale_alpha(range = c(1, 0.2)) +
        guides(shape = guide_legend(ncol = 1), color = guide_legend(override.aes=list(shape=19, size = 4))) +
        scale_colour_gradient(low = random_lightblue, high = "black") +
        theme_classic() + 
        theme(legend.key.size = unit(0.5, 'cm'), text = element_text(size = 10),
              plot.title = element_text(hjust = 0.5, size = 12), 
              panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + 
        labs(x = "GO Term", y = "Fold Enrichment", size = "Genes", color = "Adjusted\nP-Value", title = label)

    filename <- paste0("Signatures/results/figures/GREAT/", label, "_ngene_", n_genes, ".png")
    png(filename, width = 8, height = 7, res = 600, units = "in")
    print(p)
    dev.off()
}

#' Plot GREAT results with annotations
#' 
plot_GREAT_anno <- function(great, n_genes, label, anno) {

    # keep only top 20 with > n_genes
    great <- great[great$Total_Genes_Annotated > n_genes,]
    great <- great[order(great$Hyper_Fold_Enrichment, decreasing = TRUE),]
    toPlot <- great[1:20,]
    toPlot$name <- factor(toPlot$name, levels=rev(toPlot$name))

    # annotation bar
    anno$Label <- gsub("Developmental process", "Developmental\nprocess", anno$Label)
    anno <- anno[anno$ARCHE == label,]
    toPlot$anno <- anno$Label #[match(toPlot$ID, anno$GO)]

    p1 <- ggplot(toPlot, aes(x = name, y = 1, fill = anno)) +
        geom_tile(color = "black") +
        coord_flip() +
        scale_fill_manual(values = c(random_blue, "white")) +
        theme_void() +
        theme(
            axis.text.y = element_text(hjust = 1),
            axis.title.y = element_text(angle = 90),
            plot.margin=unit(c(1.03,-0.1,1.3,1), "cm")
        ) +
        labs(x = "GO Term", fill = "Related Pathways")

    # lollipop
    x <- max(toPlot$Hyper_Fold_Enrichment) + 0.25

    p2 <- ggplot(toPlot, aes(x = name, y = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH)) +
        geom_point(aes(size = Total_Genes_Annotated), shape = 19) + 
        geom_segment(aes(x = name, xend = name, y = 0, yend = Hyper_Fold_Enrichment, color = Hyper_Adjp_BH), size = 1) +
        coord_cartesian(clip = "off") + 
        coord_flip() +
        scale_y_continuous(limits = c(0, x), expand = c(0,0)) +
        guides(shape = guide_legend(ncol = 1), color = guide_legend(override.aes=list(shape=19, size = 4))) +
        scale_colour_gradient(low = random_lightblue, high = "black") +
        theme_classic() + 
        theme(
            legend.key.size = unit(0.5, 'cm'), text = element_text(size = 10),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 12), 
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            plot.margin=unit(c(0.5,1,0.5,0), "cm")
        ) + 
        labs(y = "Fold Enrichment", size = "Number of Genes", color = "Adjusted P-Value", title = paste(label, "      "))

    # extract legends
    l1 <- as_ggplot(get_legend(p1))
    l2 <- as_ggplot(get_legend(p2))
    p1 <- p1+theme(legend.position = "none")
    p2 <- p2+theme(legend.position = "none")

    filename <- paste0("Signatures/results/figures/GREAT/", label, "_annotated_ngene_", n_genes, ".png")
    png(filename, width = 10, height = 7, res = 600, units = "in")
    print(grid.arrange(p1, p2, l1, l2, ncol = 7, nrow = 6,
        layout_matrix = rbind(c(1,1,1,1,2,2,NA),
                              c(1,1,1,1,2,2,3),
                              c(1,1,1,1,2,2,3),
                              c(1,1,1,1,2,2,4),
                              c(1,1,1,1,2,2,NA),
                              c(1,1,1,1,2,2,NA))))
    dev.off()

}