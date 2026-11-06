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

#' Plot subtype enrichment of ARCHE scores
#' 
plot_subtype_enrichment <- function(scores, meta, label) {

    scores$ARCHE <- rownames(scores)
    toPlot <- reshape2::melt(scores)
    toPlot$Subtype <- meta$subtype[match(toPlot$variable, meta$sampleid)]

    p1 <- ggplot(toPlot, aes(x = ARCHE, y = value, group = variable, fill = Subtype)) +
        geom_line(aes(color = Subtype), alpha = 0.5) +
        geom_point(shape = 21, size = 2) +
        scale_fill_manual(values = subtype_pal) +
        scale_color_manual(values = subtype_pal) +
        theme_classic() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.key.size = unit(0.5, "mm")
        ) + 
        labs(y = "ARCHE Score")

    p2 <- ggplot(toPlot, aes(x = ARCHE, y = 1, fill = ARCHE)) +
        geom_tile() +
        scale_fill_manual(values = ARCHE_pal) +
        theme_void() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(size = 8)
        )

    p <- p1 / p2 + plot_layout(height = c(20, 1))
    
    filename <- paste0("data/results/figures/2.5-PreclinicalVal/samplescoring/", label, ".png")
    ggsave(filename, p, width = 5, height = 3)

}

#' Plot ARCHE scores per ARCHE ~ sample
#' 
plot_ARCHE_scores_heatmap <- function(df, label, meta, folder = "3-DataExploration") {

    # unnormalized
    sumdev <- colSums(abs(df))

    # set colours for plotting
    lim <- max(c(abs(min(df)), max(df)))
    score_pal = colorRamp2(seq(-lim, lim, length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))
    count_pal <- colorRamp2(seq(min(sumdev), max(sumdev), length = 3), c("#DFDFDF", "#989898", "#202020"))

    ha <- HeatmapAnnotation(
        Subtype = meta[match(colnames(df), meta$sampleid),]$subtype,
        Tech = meta[match(colnames(df), meta$sampleid),]$tech,
        SumDevs = sumdev,
        col = list(Subtype = subtype_pal, Tech = tech_pal, SumDevs = count_pal))

    filename <- paste0("data/results/figures/", folder, "/ARCHEheatmaps/", label, "_ARCHE_scores_unnorm.png")
    #png(filename, width = 10, height = 4, res = 600, units = "in")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(
        Heatmap(df, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
            row_names_gp = gpar(fontsize = 10), top_annotation = ha)
    )
    dev.off()

    # normalize
    df <- znorm(df)

    # set colours for plotting
    lim <- max(c(abs(min(df)), max(df)))
    score_pal = colorRamp2(seq(-lim, lim, length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))

    ha <- HeatmapAnnotation(
        Subtype = meta[match(colnames(df), meta$sampleid),]$subtype,
        Tech = meta[match(colnames(df), meta$sampleid),]$tech,
        SumDevs = sumdev,
        col = list(Subtype = subtype_pal, Tech = tech_pal, SumDevs = count_pal))

    filename <- paste0("data/results/figures/", folder, "/ARCHEheatmaps/", label, "_ARCHE_scores_norm.png")
    #png(filename, width = 10, height = 4, res = 600, units = "in")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(
        Heatmap(df, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
            row_names_gp = gpar(fontsize = 10), top_annotation = ha)
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
    png(filename, width = 5, height = 8, res = 600, units = "in")
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
    png(filename, width = 5, height = 8, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = ARCHE, y = variable, fill = value)) +
        geom_tile() +
        facet_wrap(~Sites) +
        theme_minimal() +
        scale_fill_gradient2(low = "#B05C58", mid = "#F8F1F8", high = "#077293") +
        labs(fill = "Z-Normalized\nScores", y = sample))
    dev.off()
}
