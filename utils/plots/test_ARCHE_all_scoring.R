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
        scale_shape_manual(values = sample_type_shapes) +
        scale_color_manual(values = subtype_pal) +
        theme_bw()


    p2 <- ggplot(pca_res, aes(x = PC1, y = PC2, color = type)) +
        geom_point(size = 2) +
        scale_shape_manual(values = sample_type_shapes) +
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

# helper function to plot rank waterfall plots of tumour ARCHE scores
plot_ARCHE_waterfall <- function(arche, scores) { 

    toPlot <- scores[order(scores[[arche]], decreasing = TRUE),]
    toPlot$sampleID <- factor(toPlot$sampleID, levels = toPlot$sampleID)

    # get mismatch samples
    to_label <- c()
    n <- nrow(scores[scores$ARCHE == arche,])
    top_n <- toPlot[1:n,]
    if (nrow(top_n[top_n$ARCHE != arche,]) > 0) {
        ott <- toPlot[n+1:nrow(toPlot),]
        to_label <- c(top_n$sampleID[top_n$ARCHE != arche],  ott$sampleID[ott$ARCHE == arche])

        toPlot$label <- ifelse(toPlot$sampleID %in% to_label, "*", "")
    } else {
        toPlot$label <- ""
    }

    p <- ggplot(toPlot, aes(x = sampleID, y = .data[[arche]], fill = ARCHE)) +
        geom_bar(stat = "identity", color = "black") +
        geom_text(aes(label = label), vjust = -0.1, size = 6) +
        scale_fill_manual(values = ARCHE_pal) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),
            axis.title.x = element_blank()
        )
    
    return(p)
}

# function to plot heatmap and waterfall plots of tumour ARCHE scores
plot_tumour_ARCHE_scores <- function(scores, dir, tumour_meta, label) {

    scores <- as.data.frame(t(scores[,colnames(scores) %in% tumour_meta$sampleid]))
    scores$ARCHE <- tumour_meta$arche[match(rownames(scores), tumour_meta$sampleid)]
    scores$sampleID <- rownames(scores)
    scores$rank <- as.character(tcga_meta$rank[match(scores$sampleID, gsub("\\.", "-", tcga_meta$Sample.Name))])

    # --- heatmap

    toPlot <- reshape2::melt(scores)
    toPlot$variable <- factor(toPlot$variable, levels = paste0("ARCHE", 6:1))
    toPlot <- toPlot[order(as.numeric(toPlot$rank)),]
    toPlot$sampleID <- factor(toPlot$sampleID, levels = unique(toPlot$sampleID))

    p1 <- ggplot(toPlot, aes(x = sampleID, y = "Assigned\nARCHE", fill = ARCHE)) +
        geom_tile() +
        scale_fill_manual(values = ARCHE_pal) +
        theme_void() +
        theme(legend.position = "none")

    p2 <- ggplot(toPlot, aes(x = sampleID, y = variable, fill = value)) +
        geom_tile() +
        scale_fill_gradientn("ARCHE\nExpression\nScore", colours = brewer.pal(9, "Blues")) +
        theme_void() +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
            axis.text.y = element_text(size = 9, hjust=0)
        )

    p <- p1 / p2 + plot_layout(heights = c(1, 7))

    filename <- paste0("data/results/figures/Misc/all_scoring/", dir, "/tumour_heatmap_", label, ".png")
    ggsave(filename, p, w = 9, h = 4)

    # --- individual waterfall plots

    p1 <- plot_ARCHE_waterfall("ARCHE1", scores)
    p2 <- plot_ARCHE_waterfall("ARCHE2", scores)
    p3 <- plot_ARCHE_waterfall("ARCHE3", scores)
    p4 <- plot_ARCHE_waterfall("ARCHE4", scores)
    p5 <- plot_ARCHE_waterfall("ARCHE5", scores)
    p6 <- plot_ARCHE_waterfall("ARCHE6", scores)

    p <- (p1 + p2) / (p3 + p4) / (p5 + p6) + plot_layout(guides="collect")

    filename <- paste0("data/results/figures/Misc/all_scoring/", dir, "/tumour_waterfall_", label, ".png")
    ggsave(filename, p, w = 14, h = 8)

}
