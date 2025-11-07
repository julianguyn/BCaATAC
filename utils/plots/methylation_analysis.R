#' Plot bvals per ARCHE
#' 
plot_bval_arches <- function(toPlot) {
    filename <- "data/results/figures/2-MolecularSigAnalysis/methylation/bval_arches.png"
    png(filename, width = 6, height = 4, res = 600, units = "in")
    print(
        ggplot(toPlot, aes(x = label, y = prop, fill = status)) +
            geom_bar(stat = "identity", position = "stack") +
            facet_wrap(~ ARCHE) +
            geom_text(aes(label = paste0(round(prop * 100), "%")), position = position_fill(vjust = 0.5)) +
            scale_fill_manual(values = c(random_blue, "gray", random_lightblue)) +
            theme_classic() +
            theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
            labs(x = "", y = "Proportion of CpG sites", fill = "Methylation Status")
    )
    dev.off()
}

#' Plot DMPs
#' 
plot_dmps_arches <- function(toPlot, counts) {
    filename <- "data/results/figures/2-MolecularSigAnalysis/methylation/dmps_arches.png"
    png(filename, width = 6, height = 4, res = 600, units = "in")
    print(
        ggplot() +
            geom_col(data = toPlot, aes(x = label, y = n, fill = anno)) +
            geom_text(data = counts, aes(x = label, y = n, label = n), vjust = -0.5) +
            facet_wrap(~ARCHE) +
            scale_y_continuous(limits = c(0, 117), expand = c(0, 0)) +
            scale_fill_manual("CpG Island\nMapping", values = dmp_pal) +
            theme_classic() + 
            theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
            labs(y = "Count", x = "")
    )
    dev.off()
}
