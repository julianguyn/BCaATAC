plot_ARCHE_scores_heatmap <- function(df, label, meta) {
    
    # z score normalization for each column (sample)
    means <- colMeans(df)
    stdev <- colSds(as.matrix(df))
    norm_df <- as.data.frame(sapply(1:ncol(df), function(i) { (df[[i]] - means[i]) / stdev[i] }))
    colnames(norm_df) <- colnames(df)
    rownames(norm_df) <- rownames(df)

    # set colours for plotting
    score_pal = colorRamp2(seq(min(norm_df), max(norm_df), length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))

    ha = HeatmapAnnotation(Subtype = meta[match(colnames(norm_df), meta$sample),]$subtype, 
                       col = list(Subtype = subtype_pal))

    filename <- paste0("data/results/figures/3-DataExploration/ARCHEheatmaps/", label, "_ARCHE_scores.png")
    png(filename, width = 8, height = 4, res = 600, units = "in")
    print(
        Heatmap(norm_df, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 10),
            top_annotation = ha)
    )
    dev.off()
}