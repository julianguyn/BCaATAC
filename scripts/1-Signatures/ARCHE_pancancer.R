# test ARCHE scoring across other cancer types

# load libraries
suppressPackageStartupMessages({
    library(matrixStats)
    library(data.table)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    library(patchwork)
})

set.seed(101)

source("utils/palettes.R")
source("utils/plots/ARCHE_scores_heatmap.R")

###########################################################
# Load in data
###########################################################

scores <- fread("data/rawdata/all_scoring/tcga_pancancer.Zscore.txt", data.table = FALSE)
rownames(scores) <- paste0("ARCHE", 1:6)
scores$V1 <- NULL

# load in pan-cancer metadata
meta <- read.csv("metadata/tcga_pancancer_meta.csv")

# load in BRCA metadata
tcga_meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
tcga_meta$ATAC.Seq.File.Name <- gsub("\\.", "-", tcga_meta$ATAC.Seq.File.Name)

###########################################################
# Map BRCA subtypes
###########################################################

meta$bca_subtype <- tcga_meta$Subtype[match(meta$sampleid, tcga_meta$ATAC.Seq.File.Name)]

###########################################################
# Plot heatmap
###########################################################

plot_heatmap <- function(scores, meta, label) {

    ha <- HeatmapAnnotation(
        Cancer = meta$cancer[match(colnames(scores), meta$sampleid)],
        Subtype = meta$bca_subtype[match(colnames(scores), meta$sampleid)],
        col = list(Cancer = cancer_type_pal, Subtype = subtype_pal),
        na_col = "white"
    )

    score_pal <- colorRamp2(seq(min(scores), max(scores), length = 3), c("#AD6A6C", "white", "#077293"))

    filename <- paste0("data/results/figures/1-Signatures/pancancer/", label, "_heatmap.png")
    cat("-----Saving plot to", filename, "\n")
    png(filename, width = 11, height = 4, res = 600, units = "in")
    print(
        Heatmap(scores, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 3),
            row_names_gp = gpar(fontsize = 10), top_annotation = ha)
    )
    dev.off()

}

plot_heatmap(scores, meta, "zscores")

###########################################################
# Format dataframe for boxplots
###########################################################

scores$ARCHE <- rownames(scores)
toPlot <- reshape2::melt(scores)
toPlot$Cancer <- meta$cancer[match(toPlot$variable, meta$sampleid)]
toPlot$Subtype <- meta$bca_subtype[match(toPlot$variable, meta$sampleid)]

toPlot$Label <- ifelse(toPlot$Cancer == "brca", "Breast Cancer", "Other Cancer")
toPlot$Cancer[toPlot$Cancer == "brca"] <- toPlot$Subtype[toPlot$Cancer == "brca"]


###########################################################
# Plot boxplots
###########################################################

toPlot$Cancer <- factor(toPlot$Cancer, levels = c(names(subtype_pal), names(cancer_type_pal)))
toPlot <- toPlot[-which(toPlot$Cancer %in% c("Normal", "Not Available")),]


# set x axis
x_max <- max(toPlot$value)
x_min <- min(toPlot$value)

plot_pancancer <- function(toPlot, legend_title) {
    p <- ggplot(toPlot, aes(y = Cancer, x = value, fill = Cancer)) +
        annotate("rect", ymin = -Inf, ymax = Inf, 
            xmin = 0,
            xmax = x_max,
            fill = "#CCCFD5", alpha = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_boxplot(outlier.size = 1) +
        scale_x_continuous(limits = c(x_min, x_max)) +
        scale_fill_manual(legend_title, values = c(cancer_type_pal, subtype_pal)) +
        facet_grid(
            cols = vars(ARCHE)
        ) +
        theme_bw() +
        theme(
            panel.spacing = unit(0, "lines"),
            strip.placement = "outside",
            strip.background = element_blank(),
            axis.title.y = element_blank()
        ) +
        labs(x = "Chromvar Zscores")
    return(p)
}

p1 <- plot_pancancer(toPlot[toPlot$Label == "Breast Cancer",], "BCa Subtype") +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
    )
p2 <- plot_pancancer(toPlot[toPlot$Label == "Other Cancer",], "Cancer Type") +
    theme(
        strip.text = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 9)
    )

p <- (p1 + p2) +
  plot_layout(heights = c(1, 1.5), guides = "collect") &
  theme(legend.position = "right")

filename <- "data/results/figures/1-Signatures/pancancer/pancancer_boxplots.png"
ggsave(filename, p, width = 6, height = 4)
