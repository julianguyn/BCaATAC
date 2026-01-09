#' Plot preclinical ARCHE scores
#' 
plot_ARCHE_score_preclinical <- function(toPlot) {
    png("data/results/figures/5-cfDNA/ARCHE_score_preclinical.png", width=150, height=100, units='mm', res = 600, pointsize=80)
    print(ggplot(toPlot, aes(x = ARCHE, y = Score, fill = Sample)) + 
        geom_bar(stat = "identity", position = position_dodge(), color = "black") + 
        theme_classic() + 
        geom_hline(yintercept = 0) + 
        scale_fill_manual("Sample Type", values = preclinical_pal) +
        theme(
            legend.key.size = unit(0.7, 'cm'),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
        ) + 
        labs(x = "ARCHE", y = "ARCHE Score"))
    dev.off()
}

#' Plot CICADA ARCHE scores
#' 
plot_ARCHE_score_CICADA <- function(toPlot, label) {
    filename <- paste0("data/results/figures/5-cfDNA/ARCHE_score_", label, ".png")
    png(filename, width=10, height=4, units='in', res = 600, pointsize=80)
    print(ggplot(toPlot, aes(x = ARCHE, y = mean_score, fill = label)) + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") + 
    facet_wrap(TF_group ~ .) +
    geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                    width = 0.2, position = position_dodge(width = 0.9)) +
    theme_classic() + 
    geom_hline(yintercept = 0) + 
    scale_y_continuous(limits = c(-0.05, 0.65), expand = c(0,0)) +
    scale_fill_manual("Sample Type", values = cohort_pal) +
    theme(
        legend.key.size = unit(0.7, 'cm'),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    ) + 
    labs(x = "ARCHE", y = "ARCHE Score"))
    dev.off()
}

#' Plot CICADA ARCHE scores
#' 
plot_ARCHE_TF_CICADA <- function(df, label) {

    colnames(df)[colnames(df) == "label"] <- "Sample_Type"

    filename <- paste0("data/results/figures/5-cfDNA/ARCHE_TR_corr", label, ".png")
    png(filename, width=6, height=3.5, units='in', res = 600, pointsize=80)
    print(ggplot(df, aes(x = Score, y = TF, color = Sample_Type, shape = Sample_Type)) + 
    geom_point(size = 2, alpha = 0.8) + 
    facet_wrap(ARCHE ~ ., nrow = 2, ncol = 3) +
    geom_hline(yintercept = 10, linetype = "dotted") +
    scale_color_manual(values = cohort_pal2) +
    theme_classic() + 
    theme(legend.key.size = unit(0.5, 'cm'),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)) + 
    labs(x = "ARCHE Score", y = "Tumour Fraction"))
    dev.off()
}

#' Function to coverage plots for each ARCHE
#' 
#' @param cov data.frame. Coverage at positions (+ site_name and sample name)
#' @param meta data.frame. Metadata file
#' @param group string. Site name label (10k, 20k, or 50k)
#' @param label string. Direcotry label
#' 
plot_coverage <- function(cov, meta, group, label) {
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
