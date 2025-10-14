#' Plot CICADA ARCHE scores
#' 
plot_ARCHE_score_CICADA <- function(toPlot) {
    png("cfDNA/results/figures/ARCHE_score_CICADA.png", width=10, height=4, units='in', res = 600, pointsize=80)
    print(ggplot(toPlot, aes(x = ARCHE, y = mean_score, fill = label)) + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") + 
    facet_wrap(TF_group ~ .) +
    geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                    width = 0.2, position = position_dodge(width = 0.9)) +
    theme_classic() + 
    geom_hline(yintercept = 0) + 
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
plot_ARCHE_TF_CICADA <- function(df) {

    colnames(df)[colnames(df) == "label"] <- "Sample_Type"

    png("cfDNA/results/figures/ARCHE_TR_corr.png", width=6, height=3.5, units='in', res = 600, pointsize=80)
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
