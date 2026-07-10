# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(matrixStats)

    
    library(patchwork)
    library(ggnewscale)
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/plots/drug_response_ccls.R")
source("utils/palettes.R")
source("utils/bca_drugs.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/ccl_benchmarks.R")
source("utils/plots/drug_response_pdx.R")

###########################################################
# Prepare metadata
###########################################################

# read in sample metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]

# get cell lines
c_meta <- meta[meta$type == "cell_line", ]

###########################################################
# Load in cell line data
###########################################################

# zscores
zscore_cells <- get_arche_scores(paste0("data/rawdata/all_scoring/cell_tcga.Zscore.txt"), c_meta)
normzs_cells <- znorm(zscore_cells)

# sumdevs
zscore_cells_sumdev <- get_arche_sumdevs(zscore_cells, "zscore_cells", plot = TRUE)
normzs_cells_sumdev <- znorm(zscore_cells_sumdev)

# get drug sensitivity data
load("data/procdata/CCLs/sensitivity_data.RData")

###########################################################
# Run drug repurposing pipeline
###########################################################

ubr1 <- drug_repurposing(ubr1_sen, "UBR1")
ubr2 <- drug_repurposing(ubr2_sen, "UBR2")
gray <- drug_repurposing(gray_sen, "GRAY")
gcsi <- drug_repurposing(gcsi_sen, "gCSI")
gdsc <- drug_repurposing(gdsc_sen, "GDSC2")
ccle <- drug_repurposing(ccle_sen, "CCLE")
ctrp <- drug_repurposing(ctrp_sen, "CTRP")

compiled_dr <- rbind(ubr1, ubr2, gray, gcsi, gdsc, ccle, ctrp)
compiled_dr <- compiled_dr[compiled_dr$rho < -0.3 & compiled_dr$diff > 0,]

save(compiled_dr, ubr1, ubr2, gray, gcsi, gdsc, ccle, ctrp,
    file = "data/procdata/CCLs/drug_repurposing/all_psets.RData")

table(compiled_dr$PSet)

###########################################################
# Run drug repurposing pipeline
###########################################################


arche = "ARCHE3"
sen <- gray_sen
drug1 = "AZD6244"
drug2 = "Trastuzumab"

tt <- t(zscore_cells_sumdev[arche,]) |> as.data.frame()
tt <- tt[order(tt[[arche]], decreasing = TRUE),,drop = FALSE]
cell_order <- rownames(tt)


toPlot <- t(sen[rownames(sen) %in% c(drug1, drug2), colnames(sen) %in% c(cell_order)]) |> as.data.frame()
toPlot$ARCHE <- tt[[arche]][match(rownames(toPlot), rownames(tt))]
toPlot <- toPlot[order(toPlot$ARCHE, decreasing = TRUE),]
toPlot <- toPlot[1:(nrow(toPlot)/2),]

toPlot <- data.frame(
    Sample = rep(rownames(toPlot), 2),
    AAC = c(toPlot[[drug1]], toPlot[[drug2]]),
    Drug = c(rep(drug1, nrow(toPlot)), rep(drug2, nrow(toPlot))),
    ARCHE = rep(toPlot$ARCHE, 2)
)

toPlot$Sample <- factor(toPlot$Sample, levels = unique(toPlot$Sample))
ggplot(toPlot, aes(x = ARCHE, y = AAC, fill = Drug)) +
    geom_point(shape = 21) +
    geom_smooth(method = "lm")
