# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
})

source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# label dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta$sampleid[meta$sampleid %in% dups] <- paste0(meta$sampleid[meta$sampleid %in% dups], " (", meta$tech[meta$sampleid %in% dups], ")")

c_meta <- meta[meta$type == "cell_line", c("filename", "sampleid", "subtype", "tech")]

# helper function to manually load in ARCHE scores
manual_ARCHE_scores <- function(path) {
    scores <- read.table(path)
    rownames(scores) <- paste0("ARCHE", 1:6)
    return(scores)
}

scores_so <- manual_ARCHE_scores("data/rawdata/cfDNA-ATAC/chromvar-samplesOnly/cfDNA_ATAC.Zscore.txt")
scores_as <- manual_ARCHE_scores("data/rawdata/cfDNA-ATAC/chromvar-allSamples/cfDNA_ATAC_all.Zscore.txt")

# standardize sample names of cell lines across all samples
colnames(scores_as) <- sub("^X", "", gsub("\\.(?!$)", "-", colnames(scores_as), perl = TRUE))
idx <- which(colnames(scores_as) %in% c_meta$filename)
colnames(scores_as)[idx] <- meta$sampleid[match(colnames(scores_as)[idx], meta$filename)]
colnames(scores_as) <- sub("ATAC_", "", sub("_ATAC", "", sub("_S.*", "", colnames(scores_as))))
scores_as <- scores_as[, order(colnames(scores_as))]

# standardize sample names of cfDNA samples only
colnames(scores_so) <- sub("ATAC_", "", sub("_ATAC", "", sub("_S.*", "", colnames(scores_so))))

###########################################################
# Create metadata for cfDNA samples
###########################################################

# cfDNA samples only
reps <- c("Rep1", "Rep2")
meta <- data.frame(
  SampleID = colnames(scores_so),
  Sample = sub("_.*", "", sub("RES.*", "", sub("PAR.*", "", colnames(scores_so)))),
  Model = c(rep("CCL", 12), "PDO", "PDX", "PDX", "PDO", rep("CCL", 8)),
  Res = "Parental",
  Rep = c(rep(reps, 6), "NA", reps, "NA", rep(reps, 4))
)
meta$Res[grep("RES", colnames(scores_so))] <- "Resistant"

# merge into cell metadata
to_merge <- data.frame(
    filename = meta$SampleID,
    sampleid = meta$SampleID,
    subtype = "LumB",
    tech = "sasha"
)
c_meta <- rbind(c_meta, to_merge)

###########################################################
# Plot heatmap
###########################################################

plot_ARCHE_scores_heatmap(scores_as, "allSamples", c_meta, folder = "5-cfDNA")

###########################################################
# Subset for cfDNA samples
###########################################################

scores_as <- scores_as[,meta$SampleID]

###########################################################
# Plot scores
# modify from utils/plots/cfdna.R
###########################################################

# format metadata
meta$Sample[meta$SampleID %in% c("CAMA1RESA", "CAMA1RESB")] <- "CAMA1_xeno"
meta$Rep <- paste0("-", meta$Rep)
meta$Rep[meta$Rep == "-NA"] <- ""
meta$Label <- paste0(meta$Sample, meta$Rep)
meta$Label <- factor(meta$Label, levels = names(preclinical_pal2))

plot_ARCHE_preclinical_comp <- function(scores, meta, label) {

  # format dataframe for plotting
  scores <- t(scores) 
  meta <- meta[match(rownames(scores), meta$SampleID),]
  scores <- cbind(scores, meta) |> as.data.frame()
  toPlot <- reshape2::melt(scores)
    
  filename <- paste0("data/results/figures/5-cfDNA/cfDNA-ATAC/", label, ".png")
  png(filename, width=8.5, height=5, units='in', res = 600, pointsize=80)
  print(ggplot(toPlot, aes(x = variable, y = value, fill = Label)) + 
      geom_bar(stat = "identity", position = position_dodge2(), color = "black") + 
      geom_hline(yintercept = 0) +
      facet_wrap(~Res, ncol = 1) +
      scale_fill_manual(values = preclinical_pal2) +
      theme_minimal() +
      theme(
        legend.key.size = unit(0.5, 'cm'),
        panel.border = element_rect()
      ) +
      labs(x = "ARCHE", y = "ARCHE Score"))
  dev.off()
}

plot_ARCHE_preclinical_comp(scores_as, meta, "cfDNA_allSamples")
plot_ARCHE_preclinical_comp(scores_so, meta, "cfDNA_only")

###########################################################
# Manually plot heatmap
# code from: plot_ARCHE_scores_heatmap()
###########################################################

plot_ARCHE_scores_heatmap_cfDNA <- function(df, meta, label) {

    # normalize scores
    df <- znorm(df)

    # set colours for plotting
    lim <- max(c(abs(min(df)), max(df)))
    score_pal = colorRamp2(seq(-lim, lim, length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))

    ha <- HeatmapAnnotation(
      Sample = meta[match(colnames(df), meta$SampleID),]$Label,
      Model = meta[match(colnames(df), meta$SampleID),]$Model,
      Res = meta[match(colnames(df), meta$SampleID),]$Res,
      col = list(Sample = preclinical_pal2, Model = model_pal, Res = res_pal)
    )

    filename <- paste0("data/results/figures/5-cfDNA/ARCHEheatmaps/", label, "_ARCHE_scores_norm.png")
    png(filename, width = 9, height = 4.5, res = 600, units = "in")
    print(
        Heatmap(df, cluster_rows = FALSE, name = "ARCHE\nScore", col = score_pal,
            column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
            row_names_gp = gpar(fontsize = 10), top_annotation = ha)
    )
    dev.off()
}

plot_ARCHE_scores_heatmap_cfDNA(scores_as, meta, "cfDNA_allSamples")
plot_ARCHE_scores_heatmap_cfDNA(scores_so, meta, "cfDNA_only")
