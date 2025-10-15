#' Plot distribution of ARCHE scores
#' 
plot_ARCHE_scores <- function(toPlot) {

    png("DataExploration/results/figures/ARCHE_scores.png", width = 9, height = 4, res = 600, units = "in")
    print(ggplot(toPlot, aes(x = Value)) + geom_density(fill = "#88A0A8") + 
            facet_grid(.~ARCHE) + xlim(c(-255, 255)) +
            theme_classic() + geom_vline(xintercept = 0, linetype = 'dashed', color = "gray") +
            theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
            labs(y = 'Density', x = 'ARCHE Score'))
    dev.off()
}

#' Plot overlapping BCa drugs in psets
#' 
plot_BCa_drugs_overlap <- function(toPlot) {

    png("DataExploration/results/figures/BCa_drugs_overlap.png", width=6, height=6, units='in', res = 600, pointsize=80)
    print(ggplot(toPlot, aes(x = PSet, y = Drug, fill = Present)) +
    geom_tile(color = "black") +
    scale_fill_manual(values = c("#3E517A", "white")) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
            panel.grid = element_blank(),
            legend.key.size = unit(0.5, 'cm')) +
    labs(fill = "Drug\nPresent"))
    dev.off()

}
