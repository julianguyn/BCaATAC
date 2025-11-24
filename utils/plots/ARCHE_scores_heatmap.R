#' Helper function to z-score normalization ARCHE scores for plotting
#' 
# z score normalization for each column (sample)
znorm <- function(df) {
    means <- colMeans(df)
    stdev <- colSds(as.matrix(df))
    norm_df <- as.data.frame(sapply(1:ncol(df), function(i) { (df[[i]] - means[i]) / stdev[i] }))
    colnames(norm_df) <- colnames(df)
    rownames(norm_df) <- rownames(df)
    return(norm_df)
}

#' Plot ARCHE scores per ARCHE ~ sample
#' 
plot_ARCHE_scores_heatmap <- function(df, label, meta) {

    # normalize
    df <- znorm(df)

    # set colours for plotting
    score_pal = colorRamp2(seq(min(df), max(df), length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))

    ha = HeatmapAnnotation(Subtype = meta[match(colnames(df), meta$sample),]$subtype, 
                       col = list(Subtype = subtype_pal))

    filename <- paste0("data/results/figures/3-DataExploration/ARCHEheatmaps/", label, "_ARCHE_scores.png")
    png(filename, width = 8, height = 4, res = 600, units = "in")
    print(
        Heatmap(df, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 10))
            ,
            #top_annotation = ha)
    )
    dev.off()
}

#' Compare subsetted and all sites
#' 
plot_ARCHE_scores_compare <- function(p2, p5, pa, sample) {

    # normalize
    p2_norm <- znorm(p2)
    p5_norm <- znorm(p5)
    pa_norm <- znorm(pa)

    # ----- Plot unnormalized

    # format dataframes
    p2$ARCHE <- sub("_.*", "", rownames(p2))
    p2 <- reshape2::melt(p2)
    p2$Sites <- "20k sites"

    p5$ARCHE <- sub("_.*", "", rownames(p5))
    p5 <- reshape2::melt(p5)
    p5$Sites <- "50k sites"

    pa$ARCHE <- sub("_.*", "", rownames(pa))
    pa <- reshape2::melt(pa)
    pa$Sites <- "All sites"

    # create toPlot
    toPlot <- rbind(p2, p5, pa)
    toPlot$ARCHE <- sub("ARCHE", "", toPlot$ARCHE)

    filename <- paste0("data/results/figures/3-DataExploration/ARCHEheatmaps/", sample, "_compare_unnorm.png")
    png(filename, width = 5, height = 6, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = ARCHE, y = variable, fill = value)) +
        geom_tile() +
        facet_wrap(~Sites) +
        theme_minimal() +
        scale_fill_gradient2(low = "#B05C58", mid = "#F8F1F8", high = "#077293") +
        labs(fill = "Unnormalized\nScores", y = sample))
    dev.off()

    # ----- Plot normalized (same code)

    # format dataframes
    p2_norm$ARCHE <- sub("_.*", "", rownames(p2_norm))
    p2 <- reshape2::melt(p2_norm)
    p2$Sites <- "20k sites"

    p5_norm$ARCHE <- sub("_.*", "", rownames(p5_norm))
    p5 <- reshape2::melt(p5_norm)
    p5$Sites <- "50k sites"

    pa_norm$ARCHE <- sub("_.*", "", rownames(pa_norm))
    pa <- reshape2::melt(pa_norm)
    pa$Sites <- "All sites"

    # create toPlot
    toPlot <- rbind(p2, p5, pa)
    toPlot$ARCHE <- sub("ARCHE", "", toPlot$ARCHE)

    filename <- paste0("data/results/figures/3-DataExploration/ARCHEheatmaps/", sample, "_compare_znorm.png")
    png(filename, width = 5, height = 6, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = ARCHE, y = variable, fill = value)) +
        geom_tile() +
        facet_wrap(~Sites) +
        theme_minimal() +
        scale_fill_gradient2(low = "#B05C58", mid = "#F8F1F8", high = "#077293") +
        labs(fill = "Z-Normalized\nScores", y = sample))
    dev.off()
}
