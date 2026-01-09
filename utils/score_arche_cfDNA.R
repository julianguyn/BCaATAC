#' Function to score ARCHEs from Griffin outputs
#' 
#' @param dir filepath. Directory of Griffin outputs
#' @param meta metadata file
#' 
score_arche_cfDNA <- function(dir, meta) {

    dir <- paste0("data/rawdata/cfDNA/", dir)

    # get all griffin files
    files <- list.files(
        dir,
        recursive = TRUE,
        pattern = "GC_corrected.coverage.tsv"
    )
    samples <- gsub("/.*", "", files)

    # create dataframe to store results
    res <- data.frame(matrix(nrow=0, ncol=3)) 
    colnames(res) <- c("Sample", "ARCHE", "Score")

    for (file in files) {
        df <- fread(paste0(dir, "/", file))
        df  <- data.frame(Sample = df$sample,
                          Label = df$site_name,
                          Score = 1-df$central_coverage)
        res <- rbind(res, df)
    }

    res$ARCHE <- sub("_.*", "", res$Label)
    res$Subset <- sub(".*_", "", res$Label)

    # add metadata variables
    res$label <- meta$time_id[match(res$Sample, meta$sample_id)]
    res$TF <- meta$metrics_tf[match(res$Sample, meta$sample_id)]
    res$pheno <- meta$pheno_id[match(res$Sample, meta$sample_id)]
    res$TF_group <- ifelse(res$TF > 10, "TF>10", "TF<10")

    return(res)
}

#' Function to summarize ARCHE scores
#' 
#' Compute mean and SE of ARCHE scores
#' @param scores df. Output of score_arche_cfDNA()
#' 
summary_ARCHE <- function(df, TF_group) {
    summary_df <- df %>%
        group_by(ARCHE, label) %>%
        summarise(
            mean_score = mean(Score, na.rm = TRUE),
            se_score = sd(Score, na.rm = TRUE) / sqrt(n())
        ) %>%
        ungroup()
    summary_df <- as.data.frame(summary_df)
    summary_df$TF_group <- TF_group
    return(summary_df)
}

#' Function to coverage plots for each ARCHE
#' 
#' @param cov data.frame. Coverage at positions (+ site_name and sample name)
#' @param group string. Site name label (10k, 20k, or 50k)
#' @param label string. Direcotry label
#' 
plot_cov <- function(cov, group, label) {
    toPlot <- cov[grep(group, cov$site_name), ]
    toPlot <- reshape2::melt(toPlot)
    toPlot$TF <- meta$metrics_tf[match(toPlot$sample, meta$sample_id)]
    toPlot$variable <- as.numeric(as.character(toPlot$variable))

    message(paste("Using the scales:", min(toPlot$value), max(toPlot$value)))

    p <- ggplot(toPlot, aes(x = variable, y = value, color = TF)) +
      geom_line(alpha = 0.65) +
      facet_wrap(~ site_name, scales = "fixed", nrow = 2) +
      coord_cartesian(ylim = c(min(toPlot$value), max(toPlot$value))) +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.key.size = unit(0.5, 'cm')
      ) +
      labs(
        y = "GC-corrected Fragment\nMidpoint Coverage",
        x = "Position relative to site (bp)"
      )

    png(paste0("data/results/figures/5-cfDNA/griffinplots/",label, "_", group,"_coverage.png"), width=9, height=6, units='in', res = 600, pointsize=80)
    print(p)
    dev.off()
}