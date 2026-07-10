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

# get drug response associations
load("data/results/data/4-DrugResponse/CCLs/ARCHE_CCLs_associations.RData")

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
compiled_dr <- compiled_dr[!is.na(compiled_dr$rho),]

save(compiled_dr, ubr1, ubr2, gray, gcsi, gdsc, ccle, ctrp,
    file = "data/procdata/CCLs/drug_repurposing/all_psets.RData")

table(compiled_dr$PSet)

###########################################################
# Get potential hits
###########################################################

pc_zscore_cells_sumdev <- pc_zscore_cells_sumdev[abs(pc_zscore_cells_sumdev$pc) > 0.4,]

for (arche in paste0("ARCHE", 1:6)) {
    lab1 <- paste0("Drug1_", arche)
    compiled_dr[[lab1]] <- NA
    subset_pc <- pc_zscore_cells_sumdev[pc_zscore_cells_sumdev$signature == arche,]
    compiled_dr[[lab1]] <- ifelse(compiled_dr$Drug1 %in% subset_pc$drug, 1, 0)

    lab2 <- paste0("Drug2_", arche)
    compiled_dr[[lab2]] <- NA
    subset_pc <- pc_zscore_cells_sumdev[pc_zscore_cells_sumdev$signature == arche,]
    compiled_dr[[lab2]] <- ifelse(compiled_dr$Drug2 %in% subset_pc$drug, 1, 0)
}

compiled_dr$Label <- paste(compiled_dr$Drug1, compiled_dr$Drug2, compiled_dr$PSet, sep = "_")
dr_res <- compiled_dr[, c(8:19)] #3683
rownames(dr_res) <- compiled_dr$Label

dr_res <- dr_res[rowSums(dr_res) > 0,] #2633
compiled_dr <- compiled_dr[compiled_dr$Label %in% rownames(dr_res),]

###########################################################
# Look for BCa drugs
###########################################################

# 234 res
bca_hits <- compiled_dr[compiled_dr$Drug1 %in% bca_drugs | compiled_dr$Drug2 %in% bca_drugs,]
bca_hits$pair <- paste(bca_hits$Drug1, bca_hits$Drug2, sep = "_")

for (arche in paste0("ARCHE", 1:6)) {
    bca_hits[[arche]] <- ifelse(
        bca_hits[[paste0("Drug1_", arche)]] == 1 &
        bca_hits[[paste0("Drug2_", arche)]] == 1,
        1, 0
    )
}

# 78 res
top_hits <- bca_hits[
    bca_hits$ARCHE1 == 1 | bca_hits$ARCHE2 == 1 | bca_hits$ARCHE3 == 1 |
    bca_hits$ARCHE4 == 1 | bca_hits$ARCHE5 == 1 | bca_hits$ARCHE6 == 1,
]

###########################################################
# Run drug repurposing pipeline
###########################################################

plot_all_drug_pair <- function(arche, pair, ylim = 1, subtype = FALSE) {

    drug1 <- sub("_.*", "", pair)
    drug2 <- sub(".*_", "", pair)

    p1 <- plot_drug_pair(arche, ubr1_sen, drug1, drug2, "UBR1", ylim, subtype)
    p2 <- plot_drug_pair(arche, ubr2_sen, drug1, drug2, "UBR2", ylim, subtype)
    p3 <- plot_drug_pair(arche, gcsi_sen, drug1, drug2, "gCSI", ylim, subtype)
    p4 <- plot_drug_pair(arche, gray_sen, drug1, drug2, "GRAY", ylim, subtype) + theme(legend.position = "none")
    p5 <- plot_drug_pair(arche, gdsc_sen, drug1, drug2, "GDSC2", ylim, subtype) + theme(legend.position = "none")
    p6 <- plot_drug_pair(arche, ccle_sen, drug1, drug2, "CCLE", ylim, subtype)
    p7 <- plot_drug_pair(arche, ctrp_sen, drug1, drug2, "CTRP", ylim, subtype) + theme(axis.title.y = element_blank())

    if (pair == "Docetaxel_Panobinostat") {
        p <- p4 + p7
        w = 8
    } else {
        p <- p4 + p5 + p7
        w = 11
    }
    if (subtype == TRUE) pair <- paste0(pair, "_bySubtype")
    filename <- paste("data/results/figures/4-DrugResponse/DrugRepurposing/", pair, ".png")
    ggsave(filename, p, width = w, height = 3.5)

}

plot_drug_pair <- function(arche, sen, drug1, drug2, label, ylim = 1, subtype = FALSE) {

    tt <- t(zscore_cells_sumdev[arche,]) |> as.data.frame()

    if(drug1 %in% rownames(sen) & drug2 %in% rownames(sen)) {
        toPlot <- t(sen[rownames(sen) %in% c(drug1, drug2), colnames(sen) %in% rownames(tt)]) |> as.data.frame()
        toPlot$ARCHE <- tt[[arche]][match(rownames(toPlot), rownames(tt))]
        toPlot <- toPlot[order(toPlot$ARCHE, decreasing = TRUE),]

        toPlot <- data.frame(
            Sample = rep(rownames(toPlot), 2),
            AAC = c(toPlot[[drug1]], toPlot[[drug2]]),
            Drug = c(rep(drug1, nrow(toPlot)), rep(drug2, nrow(toPlot))),
            ARCHE = rep(toPlot$ARCHE, 2)
        )
        toPlot$Subtype <- c_meta$subtype[match(toPlot$Sample, c_meta$sampleid)]

        if(subtype == FALSE) {
            p <- ggplot(toPlot, aes(x = ARCHE, y = AAC, fill = Drug)) +
                geom_smooth(method = "lm", color = "black", show.legend = FALSE) +
                geom_point(shape = 21, size = 3) +
                scale_fill_manual(values = c("#8E5572", "#77ACA2")) +
                guides(fill = guide_legend(override.aes = list(size = 3))) +
                ylim(c(0, ylim)) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5)) +
                labs(title = label, x = paste(arche, "Score"), y = "Drug Response (AAC)")
        } else {
            p <- ggplot() +
                geom_point(
                    data = toPlot, aes(x = ARCHE, y = AAC, fill = Drug),
                    shape = 22, size = 5, alpha = 0.8, color = "white") +
                scale_fill_manual(values = c("#8E5572", "#77ACA2")) +
                new_scale_fill() +
                geom_point(
                    data = toPlot,
                    aes(x = ARCHE, y = AAC, fill = Subtype),
                    shape = 21, size = 3) +
                scale_fill_manual(values = subtype_pal) +
                guides(fill = guide_legend(override.aes = list(size = 4))) +
                ylim(c(0, ylim)) +
                theme_bw() +
                theme(
                    plot.title = element_text(hjust = 0.5),
                    legend.key.size = unit(0.5, 'cm')
                ) +
                labs(title = label, x = paste(arche, "Score"), y = "Drug Response (AAC)")
        }
        return(p)
    } else {
        return(NULL)
    }
    
}

plot_all_drug_pair("ARCHE5", "Docetaxel_Panobinostat")
plot_all_drug_pair("ARCHE5", "Docetaxel_Panobinostat", subtype = TRUE)
plot_all_drug_pair("ARCHE4", "Tamoxifen_Selumetinib", subtype = TRUE, ylim = 0.4)
plot_all_drug_pair("ARCHE4", "Tamoxifen_Selumetinib", ylim = 0.4)
