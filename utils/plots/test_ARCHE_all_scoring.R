# helper function to get sum of magnitude deviations
get_devs <- function(df, dir, meta, label, lim) {

    devs <- colSums(abs(df)) |> as.data.frame()
    colnames(devs) <- "Sum"
    devs$Type <- meta$type[match(rownames(devs), meta$sampleid)]
    n <- paste0(as.character(nrow(devs[devs$Sum > lim,])), "/", nrow(devs))
    y <- switch(label, zscore = 18, deviations = 32)

    p <- ggplot(devs, aes(x = Sum, fill = Type)) +
        geom_histogram(color = "black", linewidth = 0.3) +
        scale_fill_manual(values = sample_type_pal) +
        geom_vline(xintercept = lim, linetype = "dashed", color = "gray") +
        scale_y_continuous(breaks = seq(0, y, by = 2)) +
        theme_bw() +
        ggtitle(paste0(label, ";  n >", lim, ": ", n))
    filename <- paste0("data/results/figures/Misc/all_scoring/", dir, "/sumdev_", label, ".png")
    ggsave(filename, p, w=6, h=4)

    to_keep <- rownames(devs[devs$Sum > lim,,drop=FALSE])
    df <- df[,to_keep]
    return(df)
}

# helper function to make heatmap
plot_ARCHE_scores_heatmap_counts <- function(df, dir, label, meta, znorm = FALSE, subset_dev = FALSE, lim = NULL) {

    # normalize
    if (znorm == TRUE) {
        cat("Normalizing\n")
        toPlot <- znorm(df)
        toPlot <- toPlot[, colSums(is.na(toPlot)) == 0]
        df <- df[,colnames(df) %in% colnames(toPlot)]
    } else {
        toPlot <- df
    }
    
    score_pal <- colorRamp2(seq(min(toPlot), max(toPlot), length = 3), c("#AD6A6C", "white", "#077293"))
    count_pal <- colorRamp2(seq(min(colSums(abs(df))), max(colSums(abs(df))), length = 3), c("#DFDFDF", "#989898", "#202020"))

    ha <- HeatmapAnnotation(
        Subtype = meta$subtype[match(colnames(df), meta$sampleid)],
        Type = meta$type[match(colnames(df), meta$sampleid)],
        SumDevs = colSums(abs(df)),
        col = list(Subtype = subtype_pal, Type = sample_type_pal, SumDevs = count_pal))

    filename <- paste0("data/results/figures/Misc/all_scoring/", dir, "/", label, "_ARCHE_scores.png")
    cat("-----Saving plot to", filename, "\n")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(
        Heatmap(toPlot, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 3),
            row_names_gp = gpar(fontsize = 10), top_annotation = ha)
    )
    dev.off()
}

# helper function to plot PCA plots
plot_pca <- function(pca_res, dir, meta) {

    rownames(pca_res) <- sub("_peaks.*", "", rownames(pca_res))
    pca_res <- as.data.frame(pca_res)
    pca_res$type <- meta$type[match(rownames(pca_res), meta$filename)]
    pca_res$subtype <- meta$subtype[match(rownames(pca_res), meta$filename)]

    p1 <- ggplot(pca_res, aes(x = PC1, y = PC2, color = subtype, shape = type)) +
        geom_point(size = 2) +
        scale_color_manual(values = subtype_pal) +
        theme_bw()


    p2 <- ggplot(pca_res, aes(x = PC1, y = PC2, color = type)) +
        geom_point(size = 2) +
        scale_color_manual(values = sample_type_pal) +
        theme_bw()

    p <- p1 + p2
    filename <- paste0("data/results/figures/Misc/all_scoring/", dir, "/pca.png")
    ggsave(filename, p, width = 9, height = 3.5)
}

# plot rank tile plot of ARCHE scores
plot_ARCHE_rank <- function(scores, arche, label, norm = FALSE) {

    if (norm == TRUE) scores <- znorm(scores)

    toPlot <- as.data.frame(t(scores[rownames(scores) == arche,]))
    toPlot$sampleid <- rownames(toPlot)
    toPlot$type <- meta$type[match(toPlot$sampleid, meta$sampleid)]
    toPlot$subtype <- meta$subtype[match(toPlot$sampleid, meta$sampleid)]

    toPlot <- toPlot[order(toPlot[[arche]], decreasing = TRUE),]
    toPlot$sampleid <- factor(toPlot$sampleid, levels = toPlot$sampleid)
    toPlot$type <- factor(toPlot$type, levels = names(sample_type_pal))

    p <- ggplot(toPlot, aes(x = sampleid, y = type, fill = subtype)) +
        geom_tile(color = "black") +
        scale_fill_manual(values = subtype_pal) +
        #scale_y_discrete(labels = c("Cells", "PDX", "Tumour")) +
        theme_void() +
        theme(
            legend.position = "none",
            axis.text.y = element_text(size = 8, hjust = 1, margin = margin(r = 3)),
            plot.title = element_text(size = 9, face = "bold"),
            panel.background = element_rect(fill = "#F0F0F0", color = NA)
        ) +
        ggtitle(paste(label, "-", arche))
    
    return(p)

}