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

#toPlot$Label <- ifelse(toPlot$Cancer == "TCGA-BRCA", "Breast Cancer", "Other Cancer")
#toPlot$Cancer[toPlot$Cancer == "TCGA-BRCA"] <- toPlot$Subtype[toPlot$Cancer == "TCGA-BRCA"]


###########################################################
# Plot boxplots
###########################################################

cancer_type_pal <- c(
    "TCGA-BRCA" = random_lightblue,
    # lung
    "TCGA-LUSC" = "#B4436C",
    "TCGA-LUAD" = "#D382A8",
    "TCGA-MESO" = "#FFD9DA",
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
    "TCGA-UCEC" = "#5C95FF",
    "TCGA-CESC" = "#0FA3B1",
    # endocrine
    "TCGA-PCPG" = "#C77573",
    "TCGA-THCA" = "#A23E48",
    "TCGA-ACC"  = "#856A5D",
    #cns
    "TCGA-LGG"  = "#7180AC",
    "TCGA-GBM"  = "#2B4570",
    # head and neck
    "TCGA-HNSC" = "#A8D0DB",
    # skin
    "TCGA-SKCM" = "#FF9E4A",
    "other" = "#DFDEDE"
)

toPlot$Cancer <- factor(toPlot$Cancer, levels = names(cancer_type_pal))
#toPlot <- toPlot[-which(toPlot$Cancer %in% c("Normal", "Not Available")),]


###########################################################
# Plot jitter plot
###########################################################

toPlot$BCa <- ifelse(toPlot$Cancer == "TCGA-BRCA", "BRCA", "Other")
toPlot$ARCHE <- factor(toPlot$ARCHE, levels = paste0("ARCHE", 1:6))

# format cancer type labels
names(cancer_type_pal) <- sub("TCGA-", "", names(cancer_type_pal))
toPlot$Label <- ifelse(toPlot$value > 50, sub("TCGA-", "", as.character(toPlot$Cancer)), "other")
toPlot$Label <- factor(toPlot$Label, levels = names(cancer_type_pal))
toPlot$Label <- droplevels(toPlot$Label)

p <- ggplot(toPlot, aes(x = ARCHE, y = value, color = Label, shape = BCa)) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "gray") +
    geom_jitter(size = 3, alpha = 0.8) +
    scale_color_manual("TCGA Cohort", values = cancer_type_pal) +
    scale_shape_manual(values = c(4, 19)) +
    theme_bw() +
    theme(
        legend.key.size = unit(0.5, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10)
    ) +
    labs(y = "Chromvar Zscores")

filename <- "data/results/figures/1-Signatures/pancancer/pancancer_jitter.png"
ggsave(filename, p, width = 11, height = 4)
