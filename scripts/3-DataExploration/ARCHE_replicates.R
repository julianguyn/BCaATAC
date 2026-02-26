# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(data.table)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
})

source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/palettes.R")
source("utils/get_data.R")
source("utils/mappings.R")

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
p_meta <- meta[meta$type == "PDX", c("filename", "sampleid", "subtype", "tech")]

# helper function to manually load in cell line replicate ARCHE scores
manual_ARCHE_scores <- function(path) {
    scores <- read.table(path)
    rownames(scores) <- paste0("ARCHE", 1:6)
    return(scores)
}

ccl <- manual_ARCHE_scores("data/rawdata/cfDNA-ATAC/chromvar-allSamples/cfDNA_ATAC_all.Zscore.txt") |> znorm()

# load in PDX arche scores
pdxs_20k <- get_arche_scores("pdxs", "k20", p_meta) |> znorm()
pdxs_50k <- get_arche_scores("pdxs", "k50", p_meta) |> znorm()
pdxs_all <- get_arche_scores("pdxs", "all", p_meta) |> znorm()

###########################################################
# Clean ccl labels
###########################################################

# standardize sample names of cell lines across all samples
colnames(ccl) <- sub("^X", "", gsub("\\.(?!$)", "-", colnames(ccl), perl = TRUE))
idx <- which(colnames(ccl) %in% c_meta$filename)
colnames(ccl)[idx] <- c_meta$sampleid[match(colnames(ccl)[idx], c_meta$filename)]
colnames(ccl) <- sub("ATAC_", "", sub("_ATAC", "", sub("_S.*", "", colnames(ccl))))
ccl <- ccl[, order(colnames(ccl))]

###########################################################
# Keep replicates
###########################################################

pdx <- pdxs_50k

# identify replicates
reps_ccl <- colnames(ccl)[grep("Rep| \\(", colnames(ccl))]
reps_pdx <- colnames(pdx)[grep("\\(", colnames(pdx))]

# keep only replicates
ccl <- ccl[,reps_ccl]
pdx <- pdx[,reps_pdx]

###########################################################
# Merge dataframes
###########################################################

# make labels
anno <- data.frame(
    sampleid = paste(
        sub("_Rep.*", "", sub(" \\(.*", "", c(colnames(ccl), colnames(pdx)))),
        rep(c("Rep1", "Rep2"), length(c(colnames(ccl), colnames(pdx)))/2),
        sep = "_"
    ),
    match = c(colnames(ccl), colnames(pdx)),
    Model = c(rep("CCL", length(colnames(ccl))), rep("PDX", length(colnames(pdx)))),
    sample = sub("_Rep.*", "", sub(" \\(.*", "", c(colnames(ccl), colnames(pdx))))
)
anno$rep_source <- ifelse(
    grepl("Rep", anno$match) == TRUE,
    "Sasha",
    "Lupien"
)

# manually fix ccl sampleids
anno$sampleid <- sub("CAMA1_PAR", "CAMA-1", anno$sampleid)
anno$sampleid <- sub("T47D_PAR", "T-47D", anno$sampleid)
anno$sampleid <- sub("ZR751_PAR", "ZR-75-1", anno$sampleid)
anno <- anno[-which(grepl("RES", anno$sampleid)),]


# merge dataframes
toPlot <- cbind(ccl, pdx) |> as.data.frame()
toPlot <- toPlot[, -which(grepl("RES", colnames(toPlot)))]
colnames(toPlot) <- anno$sampleid[match(colnames(toPlot), anno$match)]

###########################################################
# Add subtype and replicate type
###########################################################

# add tech and subtype
anno$tech <- meta$tech[match(anno$match, meta$sampleid)]
anno$Subtype <- meta$subtype[match(anno$match, meta$sampleid)]
anno$tech[is.na(anno$tech)] <- "sasha"
anno$Subtype[is.na(anno$Subtype)] <- "LumB"

# add replicate type
anno <- anno %>%
  group_by(sample) %>%
  mutate(
    Replicate = case_when(
      tech == "sasha"                          ~ "technical replicate",
      any(tech == "nergiz") & tech == "nergiz" ~ "single-end replicate",
      any(tech == "nergiz") & tech != "nergiz" ~ "paired-end replicate",
      any(tech == "tina") & any(tech == "komal") ~ "biological replicate"
    )
  ) %>%
  ungroup () %>%
  mutate(
    rep_order_temp = case_when(
        Replicate == "biological replicate" ~ 1,
        Replicate == "single-end replicate" ~ 2,
        Replicate == "paired-end replicate" ~ 2,
        Replicate == "technical replicate" ~ 3
    )
  )

# remove samples
rm <- c("BT-549", "HCC1395", "HCC1806", "HCC70", "Hs 578T", "MDA-MB-436", "REF003", "69196", "46962", "66684", "MDA-MB-231", "48602", "REF001", "T47D_PAR")
anno <- anno[-which(anno$sample %in% rm),]

# order
anno <- anno %>%
  mutate(
    rep_order = ifelse(tech == "nergiz", 1, 2),
    type_order = ifelse(Model == "PDX", 1, 2),
  ) %>%
  arrange(type_order, rep_order_temp, sample, rep_order) %>%
  dplyr::select(-type_order, -rep_order)
anno$sampleid <- paste0(
    sub("Rep.*", "", anno$sampleid),
    rep(c("Rep1", "Rep2"), length(anno$sampleid)/2)
)
anno$sampleid <- factor(anno$sampleid, levels = anno$sampleid)

anno_toPlot <- anno %>%
    dplyr::select(sampleid, Model, Subtype, Replicate) %>%
    pivot_longer(
        cols =-sampleid,
        names_to = "Variable",
        values_to = "Value"
    )
anno_toPlot$Variable <- factor(anno_toPlot$Variable, levels = c("Replicate", "Subtype", "Model"))
anno_toPlot$Value <- factor(
    anno_toPlot$Value,
    c(names(model_pal), names(replicate_pal), "LumB", "Basal", "TNBC", "TPBC")
)

###########################################################
# Plot ARCHE scores
###########################################################

# format dataframe for plotting
toPlot <- t(toPlot) |> as.data.frame()
toPlot$sampleid <- rownames(toPlot)
toPlot <- reshape2::melt(toPlot)
toPlot <- toPlot[toPlot$sampleid %in% as.character(anno$sampleid),]
toPlot$sampleid <- factor(toPlot$sampleid, levels = anno$sampleid)
toPlot$y_label <- ifelse(
    grepl("Rep1", toPlot$sampleid),
    "",
    sub("_Rep2", "", toPlot$sampleid)
)

# line positions between replicates
hline_positions <- seq(2.5, nlevels(toPlot$sampleid) - 0.5, by = 2)

p1 <- ggplot(toPlot, aes(x = variable, y = sampleid, fill = value)) + 
    geom_tile() + 
    geom_hline(yintercept = hline_positions, linewidth = 0.5) +
    scale_y_discrete(labels = setNames(toPlot$y_label, toPlot$sampleid)) +
    scale_fill_gradient2(low = "#C3BFCC", mid = "#F8F1F8", high = "#077293") +
    theme_minimal() +
    theme(
        legend.key.size = unit(0.5, 'cm'),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
    ) +
    labs(fill = "ARCHE\nScore")

p2 <- ggplot(anno_toPlot, aes(x = Variable, y = sampleid, fill = Value)) +
    geom_tile() +
    geom_hline(yintercept = hline_positions, linewidth = 0.5) +
    scale_fill_manual(
        "Subtype\nModel\nReplicate",
        values = c(subtype_pal, model_pal, replicate_pal)
    ) +
    theme_minimal() +
    theme(
        legend.key.size = unit(0.5, 'cm'),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
    )

filename <- paste0("data/results/figures/3-DataExploration/replicates_p1.png")
png(filename, width=4, height=6, units='in', res = 600, pointsize=80)
p1
dev.off()

filename <- paste0("data/results/figures/3-DataExploration/replicates_p2.png")
png(filename, width=2.5, height=6, units='in', res = 600, pointsize=80)
p2
dev.off()
