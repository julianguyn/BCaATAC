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
meta <- fread("metadata/GDC_identifiers_all.tsv", data.table = FALSE)
meta$sample_match <- sub("_atac.*", "", meta$'File Name')
meta <- meta[,c('sample_match', 'Project')] |> unique()
colnames(meta) <- c("sampleid", "cancer")

# load in BRCA metadata
tcga_meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
tcga_meta$ATAC.Seq.File.Name <- gsub("\\.", "-", tcga_meta$ATAC.Seq.File.Name)

###########################################################
# Map BRCA subtypes
###########################################################

meta$bca_subtype <- tcga_meta$Subtype[match(meta$sampleid, tcga_meta$ATAC.Seq.File.Name)]

###########################################################
# Format dataframe for boxplots
###########################################################

scores$ARCHE <- rownames(scores)
toPlot <- reshape2::melt(scores)
toPlot$Cancer <- meta$cancer[match(toPlot$variable, meta$sampleid)]
toPlot$Subtype <- meta$bca_subtype[match(toPlot$variable, meta$sampleid)]

toPlot$Label <- ifelse(toPlot$Cancer == "TCGA-BRCA", "Breast Cancer", "Other Cancer")
toPlot$Cancer[toPlot$Cancer == "TCGA-BRCA"] <- toPlot$Subtype[toPlot$Cancer == "TCGA-BRCA"]


###########################################################
# Plot boxplots
###########################################################

cancer_type_pal <- c(
    # lung
    "TCGA-LUSC" = "#7D8491",
    "TCGA-LUAD" = "#C0C5C1",
    "TCGA-MESO" = "#BCBDB8",
    # gastointestinal
    "TCGA-STAD" = "#412234",
    "TCGA-COAD" = "#6D466B",
    "TCGA-ESCA" = "#B49FCC",
    "TCGA-LIHC" = "#73689B",
    "TCGA-CHOL" = "#6E5885",
    # Genitourinary
    "TCGA-KIRP" = "#353D2F", 
    "TCGA-KIRC" = "#515B3A",
    "TCGA-BLCA" = "#607D56",
    "TCGA-PRAD" = "#6BA368",
    # soft tissue
    "TCGA-TGCT" = "#F2DC5D",
    #gynecologic
    "TCGA-UCEC" = "#60463B",
    "TCGA-CESC" = "#856A5D",
    # endocrine
    "TCGA-PCPG" = "#C77573",
    "TCGA-THCA" = "#A23E48",
    "TCGA-ACC"  = "#A43F34",
    #cns
    "TCGA-LGG"  = "#7180AC",
    "TCGA-GBM"  = "#2B4570",
    # head and neck
    "TCGA-HNSC" = "#A8D0DB",
    #"TCGA-BRCA" = "white",
    # skin
    "TCGA-SKCM" = "#A37A74"
)

toPlot$Cancer <- factor(toPlot$Cancer, levels = c(names(subtype_pal), names(cancer_type_pal)))
toPlot <- toPlot[-which(toPlot$Cancer %in% c("Normal", "Not Available")),]


# set x axis
x_max <- max(toPlot$value)
x_min <- min(toPlot$value) - 40

plot_pancancer <- function(toPlot, legend_title) {
    p <- ggplot(toPlot, aes(y = Cancer, x = value, fill = Cancer)) +
        annotate("rect", ymin = -Inf, ymax = Inf, 
            xmin = 0,
            xmax = x_max,
            fill = "#dddfe2", alpha = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
        geom_boxplot(outlier.size = 0.7, linewidth = 0.3, width = 0.5) +
        scale_x_continuous(limits = c(x_min, x_max)) +
        scale_fill_manual(legend_title, values = c(cancer_type_pal, subtype_pal)) +
        facet_grid(
            cols = vars(ARCHE)
        ) +
        theme_classic() +
        theme(
            panel.spacing = unit(0, "lines"),
            strip.placement = "outside",
            strip.background = element_blank(),
            axis.text.y = element_text(size = 7),
            axis.title.y = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
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
  plot_layout(heights = c(1, 6), guides = "collect") &
  theme(legend.position = "right")

filename <- "data/results/figures/1-Signatures/pancancer/pancancer_boxplots.png"
ggsave(filename, p, width = 7, height = 7)
