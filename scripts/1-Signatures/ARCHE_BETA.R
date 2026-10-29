# BETA results

# load libraries
suppressPackageStartupMessages({
  library(patchwork)
})

source("utils/palettes.R")

###########################################################
# Load in uptargets
###########################################################

dir <- "data/results/data/1-Signatures/BETA/"

a2 <- read.table(paste0(dir, "ARCHE2", "_uptarget.txt"), header = TRUE)
a3 <- read.table(paste0(dir, "ARCHE3", "_uptarget.txt"), header = TRUE)
a4 <- read.table(paste0(dir, "ARCHE4", "_uptarget.txt"), header = TRUE)
a5 <- read.table(paste0(dir, "ARCHE5", "_uptarget.txt"), header = TRUE)
a6 <- read.table(paste0(dir, "ARCHE6", "_uptarget.txt"), header = TRUE)


###########################################################
# Load in shortlist from deg_analysis.R (TODO: put this into a script)
###########################################################

shortlist <- list(
  ARCHE1 = c("NPY1R", "RERG", "CITED1"),
  ARCHE2 = c("FOXC1", "ELF5", "VGLL1", "SOX10", "KRT5", "KRT17", "CRYAB", "ACTR3B"),
  ARCHE3 = c("ERBB2", "C17orf37", "GSDMB", "PPP1R1B", "S100A8", "DMBT1"),
  ARCHE4 = c("IGF1R", "TRH"),
  ARCHE5 = c("SPIB", "CD274", "LAG3", "IDO1", "GZMB", "IFNG", "CXCL10", "VCAM1"),
  ARCHE6 = c("FOXJ1", "CDC20B", "IGFALS", "CXCL14")
)

###########################################################
# BETA plot
###########################################################

plot_BETA <- function(beta, arche) {
    genes <- unlist(unname(shortlist[arche]))
    toPlot <- beta[beta$padj < 0.05,]
    toPlot$label <- NA_character_

    # get shortlist genes
    toPlot$label[toPlot$GeneSymbol %in% genes & toPlot$log2FC >= 1] <- "DEG-Function"

    # get top 10 RP genes
    top_rp_genes <- head(
        toPlot$GeneSymbol[order(toPlot$reg_potential, decreasing = TRUE)],
        10
    )

    # Assign top RP and overlap labels
    toPlot$label[toPlot$GeneSymbol %in% top_rp_genes &
                is.na(toPlot$label) & toPlot$log2FC >= 1] <- "Top RP"

    toPlot$label[toPlot$GeneSymbol %in% top_rp_genes &
                toPlot$GeneSymbol %in% genes &
                toPlot$log2FC >= 1] <- "Overlap"

    p <- ggplot(toPlot, aes(x = reg_potential, y = log2FC)) +
        geom_point(aes(color = label), size = 0.8) +
        geom_hline(yintercept = 1, linetype = "dashed",
           color = "gray40", linewidth = 0.4) +
        geom_point(
            data = subset(toPlot, is.na(label)),
            color = "gray", size = 0.6
        ) +
        geom_label_repel(
            data = subset(toPlot, !is.na(label)),
            aes(label = GeneSymbol, color = label),
            show.legend = FALSE,
            size = 3,
            fill = alpha("white", 0.7),
            linewidth = 0,
            label.padding = unit(0.1, "lines"),
            max.overlaps = Inf,
            min.segment.length = 0,
            segment.size = 0.3,
            box.padding = 0.4
        ) +
        scale_color_manual(values = c(
            "DEG-Function" = "#5E4C5A",
            "Top RP" = "#55917F",
            "Overlap" = "#D99B42"
        ), na.translate = FALSE) +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 11),
            legend.key.size = unit(1, "mm"),
            axis.title.x = element_text(size = 9),
            axis.title.y = element_text(size = 9)
        ) +
        labs(
            x = "Regulatory Potential Score",
            y = "log2 Fold Change",
            title = paste(arche),
            color = NULL
        )

    return(p)
}

###########################################################
# Plot BETA results
###########################################################

p2 <- plot_BETA(a2, "ARCHE2")
p3 <- plot_BETA(a3, "ARCHE3")
p4 <- plot_BETA(a4, "ARCHE4") + theme(legend.position = "none")
p5 <- plot_BETA(a5, "ARCHE5") + theme(legend.position = "none")
p6 <- plot_BETA(a6, "ARCHE6") + theme(legend.position = "none")


p <- p2 / p3 / p4 / p5 / p6 + plot_layout(ncol = 1, guides = "collect")
ggsave("data/results/figures/1-Signatures/BETA/scatterplots.png", p, width = 4, height = 12)
