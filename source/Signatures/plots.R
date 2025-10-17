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
