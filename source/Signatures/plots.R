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
    png("Signatures/results/figures/ARCHE_peakInfo.png", width = 5, height = 5, res = 600, units = "in")
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
plot_ATAC_Upset <- function(m) {

    png("Signatures/results/figures/ATAC_Upset.png", width = 11, height = 5, res = 600, units = "in")
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
            scale_fill_manual(values = brewer.pal(11, name = "Paired")) +
            theme_minimal() + 
            labs(y = "Percentage (%)")
        )
    dev.off()
}
