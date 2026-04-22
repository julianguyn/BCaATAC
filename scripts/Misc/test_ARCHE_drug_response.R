# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(ggplot2)
    library(ggh4x)
    library(reshape2)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(dplyr)
    library(readxl)
    library(data.table)
    library(patchwork)
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/compute_drug_response.R")
source("utils/plots/drug_response_ccls.R")
source("utils/palettes.R")
source("utils/bca_drugs.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/plots/drug_response_pdx.R")

###########################################################
# Prepare metadata
###########################################################

# read in cell metadata
meta <- read.csv("metadata/lupien_metadata.csv")

# remove nergiz dups
dups <- meta$sampleid[duplicated(meta$sampleid)]
meta <- meta[!(meta$sampleid %in% dups & meta$tech == "nergiz"), ]

# remove second dups
dups <- meta$sampleid[duplicated(meta$sampleid)] # only PDXs have dups
meta_t <- meta[!(meta$sampleid %in% dups & meta$tech == "tina"), ]
meta_k <- meta[!(meta$sampleid %in% dups & meta$tech == "komal"), ]

###########################################################
# Load in cell line data
###########################################################

c_meta <- meta[meta$type == "cell_line", ]

# load in arche zscores and normalize
zscore_cells <- get_arche_scores("cells", "k50", c_meta)
normzs_cells <- znorm(zscore_cells)

# load in arche deviations and normalize
deviat_cells <- get_arche_devs("cells", "50k", c_meta)
normdv_cells <- znorm(deviat_cells)

# get drug sensitivity data
load("data/procdata/CCLs/sensitivity_data.RData")

###########################################################
# Try removing low deviation samples
###########################################################

# helper function to get sum of magnitude deviations
get_devs <- function(df, label, lim1, lim2) {

    devs <- colSums(abs(df)) |> as.data.frame()
    colnames(devs) <- "Sum"
    n1 <- paste0(as.character(length(devs[devs$Sum > lim1,])), "/", nrow(devs))
    n2 <- paste0(as.character(length(devs[devs$Sum > lim2,])), "/", nrow(devs))

    p <- ggplot(devs, aes(x = Sum)) +
        geom_histogram(fill = random_lightblue, color = "black", linewidth = 0.3) +
        geom_vline(xintercept = c(lim1, lim2), linetype = "dashed", color = "gray") +
        geom_text(stat = "bin",
            aes(label = after_stat(count)),
            vjust = -0.5, size = 3) +
        theme_bw() +
        ggtitle(paste0(label, ";  n >", lim1, ":", n1, ";  n >", lim2, ":", n2))
    filename <- paste0("data/results/figures/Misc/sumdevs/", label, ".png")
    ggsave(filename, p, w=5, h=4)

    to_keep <- rownames(devs[devs$Sum > lim1,,drop=FALSE])
    df <- df[,to_keep]
    return(df)
}

# get samples above threshold
zscore_cells_sumdev <- get_devs(zscore_cells, "zscore_cells", 150, 200)
deviat_cells_sumdev <- get_devs(deviat_cells, "deviat_cells", 0.5, 0.6)

# normalize
normzs_cells_sumdev <- znorm(zscore_cells_sumdev)
normdv_cells_sumdev <- znorm(deviat_cells_sumdev)

###########################################################
# Compute PC of ARCHE-drug associations in cell lines
###########################################################

# helper function to compute arche associations across all psets
arche_pc <- function(scores, label) {
    ubr1_PC <- compute_pc(scores, ubr1_sen, "UBR1")
    ubr2_PC <- compute_pc(scores, ubr2_sen, "UBR2")
    gray_PC <- compute_pc(scores, gray_sen, "GRAY")
    gcsi_PC <- compute_pc(scores, gcsi_sen, "gCSI")
    gdsc_PC <- compute_pc(scores, gdsc_sen, "GDSC2")
    ctrp_PC <- compute_pc(scores, ctrp_sen, "CTRP")
    ccle_PC <- compute_pc(scores, ccle_sen, "CCLE")

    # compile results
    PC_res <- rbind(ubr1_PC, ubr2_PC, gray_PC, gcsi_PC, gdsc_PC, ctrp_PC, ccle_PC)
    PC_res$Label <- label
    return(PC_res)
}

# zscores
pc_zscore_cells <- arche_pc(zscore_cells, "zscore")
pc_normzs_cells <- arche_pc(normzs_cells, "normzscr")

# deviations
pc_deviat_cells <- arche_pc(deviat_cells, "deviat")
pc_normdv_cells <- arche_pc(normdv_cells, "normdev")

# subsetted zscores
pc_zscore_cells_sumdev <- arche_pc(zscore_cells_sumdev, "zscore_sumdev")
pc_normzs_cells_sumdev <- arche_pc(normzs_cells_sumdev, "normzscr_sumdev")

# subsetted deviations
pc_deviat_cells_sumdev <- arche_pc(deviat_cells_sumdev, "deviat_sumdev")
pc_normdv_cells_sumdev <- arche_pc(normdv_cells_sumdev, "normdev_sumdev")

###########################################################
# Compile cell line results
###########################################################

cell_toPlot <- rbind(
    pc_zscore_cells, pc_normzs_cells,
    pc_deviat_cells, pc_normdv_cells,
    pc_zscore_cells_sumdev, pc_normzs_cells_sumdev,
    pc_deviat_cells_sumdev, pc_normdv_cells_sumdev
)
cell_toPlot$drug[cell_toPlot$drug == "945"] <- "CFI-400945"
cell_toPlot$pairs <- sub("_945", "_CFI-400945", cell_toPlot$pairs)

###########################################################
# Plot cell results
###########################################################

# plot cells
plot_cells <- function(pair) {
    subset <- cell_toPlot[cell_toPlot$pairs == pair,]
    subset$sig <- ifelse(
        subset$FDR < 0.1,
        ifelse(subset$FDR < 0.05, 'FDR < 0.05', 'FDR < 0.1'),
        'FDR >= 0.1')
    subset$text <- ifelse(
        subset$FDR < 0.1, subset$pc, NA
    )

    subset$Label <- factor(subset$Label, levels = c(
        "zscore", "deviat", "normzscr", "normdev",
        "zscore_sumdev", "deviat_sumdev",
        "normzscr_sumdev", "normdev_sumdev"
        ))

    p <- ggplot(subset, aes(x = pset, y = Label, fill = pc, shape = sig, size = -log(FDR))) +
        geom_point() +
        geom_text(aes(label = round(text, 2)), color = "black", size = 2.5) +
        scale_shape_manual(values = c(21, 22, 23)) +
        scale_size(range = c(2, 12)) +
        scale_fill_gradient2(
            low = "#BC4749",
            high = "#689CB0",
            mid = "#C2BBC9",
            limits = c(-0.8, 0.8)
        ) +
        theme_bw() +
        theme(legend.key.size = unit(0.3, 'cm'),
            legend.title = element_text(size = 10)) +
        ggtitle(pair)

    pair <- sub(" .*", "", pair)
    filename <- paste0("data/results/figures/Misc/ARCHE_dr_test/cells_", pair, ".png")
    ggsave(filename, p, w = 6, h = 4)
}

# ----------------- ARCHE1 -----------------
# ABT-263: Bcl-2 inhibitor
plot_cells("ARCHE1_navitoclax:pluripotin (1:1 mol/mol)")
plot_cells("ARCHE1_JQ-1:navitoclax (2:1 mol/mol)")
plot_cells("ARCHE1_navitoclax:PLX-4032 (1:1 mol/mol)")
plot_cells("ARCHE1_BRD-K27188169:navitoclax (2:1 mol/mol)")
plot_cells("ARCHE1_navitoclax:birinapant (1:1 mol/mol)")
plot_cells("ARCHE1_selumetinib:navitoclax (8:1 mol/mol)")
plot_cells("ARCHE1_AZD5991") # Mcl-1 inhibitor

# AZD-5305: PARPi
plot_cells("ARCHE1_Olaparib")

# CDX-011 & Eribulin: MT inhibitors
plot_cells("ARCHE1_Paclitaxel")
plot_cells("ARCHE1_Ixabepilone")

# Everolimus: mTORi
plot_cells("ARCHE1_Sirolimus")
plot_cells("ARCHE1_Temsirolimus")

# SG: TROP2 target, payload = TOPIi
plot_cells("ARCHE1_Teniposide") #TOPIIi

# PI3K inhibitors (no PDX match)
plot_cells("ARCHE1_Alpelisib")
plot_cells("ARCHE1_BYL719")
plot_cells("ARCHE1_Taselisib")
plot_cells("ARCHE1_Pictilisib")

# AKT inhibitors (no PDX match)
plot_cells("ARCHE1_Ipatasertib")
plot_cells("ARCHE1_Afuresertib")
plot_cells("ARCHE1_Uprosertib")
plot_cells("ARCHE1_MK-2206")

# MEK inhibitors (no PDX match)
plot_cells("ARCHE1_Trametinib")
plot_cells("ARCHE1_Selumetinib")

# Platinum analogues (no PDX match)
plot_cells("ARCHE1_Carboplatinum")
plot_cells("ARCHE1_Oxaliplatin")

# CDK4/6 inhibitors (no PDX match)
plot_cells("ARCHE1_Ribociclib")


# ----------------- ARCHE2 -----------------
# ABT-263: Bcl-2 inhibitor
plot_cells("ARCHE2_navitoclax:gemcitabine (1:1 mol/mol)")
plot_cells("ARCHE2_alisertib:navitoclax (2:1 mol/mol)")
plot_cells("ARCHE2_doxorubicin:navitoclax (2:1 mol/mol)")
plot_cells("ARCHE2_decitabine:navitoclax (2:1 mol/mol)")
plot_cells("ARCHE2_navitoclax:pluripotin (1:1 mol/mol)")
plot_cells("ARCHE2_selumetinib:navitoclax (8:1 mol/mol)")
plot_cells("ARCHE2_navitoclax:MST-312 (1:1 mol/mol)")
plot_cells("ARCHE2_navitoclax:piperlongumine (1:1 mol/mol)")
plot_cells("ARCHE2_navitoclax:PLX-4032 (1:1 mol/mol)")
plot_cells("ARCHE2_AZD5991") # Mcl-1 inhibitor

# Everolimus: mTORi
plot_cells("ARCHE2_Temsirolimus")
plot_cells("ARCHE2_Ku-0063794")

# Fluvastatin: HMG-CoA reductase inhibitor / statin
plot_cells("ARCHE2_Lovastatin")
plot_cells("ARCHE2_Simvastatin")

# Topotecan (TOPIi)
plot_cells("ARCHE2_SN-38")
plot_cells("ARCHE2_Irinotecan")
plot_cells("ARCHE2_TPT")
plot_cells("ARCHE2_Topotecan")

# Chemotherapies (no PDX match)
plot_cells("ARCHE2_Paclitaxel")
plot_cells("ARCHE2_Docetaxel")
plot_cells("ARCHE2_Doxorubicin")
plot_cells("ARCHE2_Epirubicin")
plot_cells("ARCHE2_Gemcitabine")
plot_cells("ARCHE2_Vinorelbine")

# Platinum (no PDX match)
plot_cells("ARCHE2_Carboplatinum")
plot_cells("ARCHE2_Cisplatin")

# PI3Ki (no PDX match)
plot_cells("ARCHE2_Alpelisib")
plot_cells("ARCHE2_BYL719")

# AKTi (no PDX match)
plot_cells("ARCHE2_Ipatasertib")
plot_cells("ARCHE2_Afuresertib")
plot_cells("ARCHE2_Uprosertib")
plot_cells("ARCHE2_AZD5363")
plot_cells("ARCHE2_MK-2206")

# MEKi (no PDX match)
plot_cells("ARCHE2_Selumetinib")

# PARPi (no PDX match)
plot_cells("ARCHE2_Olaparib")

# SERDs (no PDX match)
plot_cells("ARCHE2_Fulvestrant")
plot_cells("ARCHE2_GDC0810")

# Lapatinib for some reason (no PDX match)
plot_cells("ARCHE2_Lapatinib")

# Under investigation (no PDX match)
plot_cells("ARCHE2_Adavosertib")
plot_cells("ARCHE2_AZD7762")


# ----------------- ARCHE3 -----------------
# ABT-263: Bcl-2 inhibitor
plot_cells("ARCHE3_WEHI-539")
plot_cells("ARCHE3_navitoclax:birinapant (1:1 mol/mol)")
plot_cells("ARCHE3_navitoclax:pluripotin (1:1 mol/mol)")
plot_cells("ARCHE3_selumetinib:navitoclax (8:1 mol/mol)")
plot_cells("ARCHE3_navitoclax:PLX-4032 (1:1 mol/mol)")
plot_cells("ARCHE3_AZD5991") # Mcl-1 inhibitor

# Paclitaxel: MT inhibitors
plot_cells("ARCHE3_Paclitaxel")
plot_cells("ARCHE3_Docetaxel")
plot_cells("ARCHE3_Ixabepilone")
plot_cells("ARCHE3_Vinorelbine")
plot_cells("ARCHE2_Mebendazole")

# Everolimus: mTORi
plot_cells("ARCHE3_Sirolimus")
plot_cells("ARCHE3_Temsirolimus")
plot_cells("ARCHE3_Ku-0063794")

# Fluvastatin: HMG-CoA reductase inhibitor / statin
plot_cells("ARCHE3_Fluvastatin")
plot_cells("ARCHE3_Lovastatin")
plot_cells("ARCHE3_Simvastatin")

# Topotecan (TOPIi)
plot_cells("ARCHE3_Irinotecan")
plot_cells("ARCHE3_TPT")

# HER2 agents (no PDX match)
plot_cells("ARCHE3_Lapatinib")
plot_cells("ARCHE3_Trastuzumab")

# Chemotherapies (no PDX match)
plot_cells("ARCHE3_Epirubicin")
plot_cells("ARCHE3_Cyclophosphamide")
plot_cells("ARCHE3_Gemcitabine")

# Platinum (no PDX match)
plot_cells("ARCHE3_Carboplatinum")
plot_cells("ARCHE3_Cisplatin")

# PI3Ki (no PDX match)
plot_cells("ARCHE3_Alpelisib")
plot_cells("ARCHE3_BYL719")
plot_cells("ARCHE3_Pictilisib")
plot_cells("ARCHE3_Taselisib")
plot_cells("ARCHE3_AZD6482")

# AKTi (no PDX match)
plot_cells("ARCHE3_Ipatasertib")
plot_cells("ARCHE3_Afuresertib")
plot_cells("ARCHE3_Uprosertib")
plot_cells("ARCHE3_AZD5363")

# MEKi (no PDX match)
plot_cells("ARCHE3_Selumetinib")
plot_cells("ARCHE3_Trametinib")
plot_cells("ARCHE3_AZD6244")

# SERDs (no PDX match)
plot_cells("ARCHE3_Fulvestrant")
plot_cells("ARCHE3_GDC0810")

# HDACi (no PDX match)
plot_cells("ARCHE3_Panobinostat")
plot_cells("ARCHE3_Vorinostat")

# Under investigation (no PDX match)
plot_cells("ARCHE3_Adavosertib")
plot_cells("ARCHE3_AZD7762")
plot_cells("ARCHE3_AZD6738")
plot_cells("ARCHE3_KU-55933")


# ----------------- ARCHE4 -----------------


# todo: individual CFI-400945 plots

######################################################################################################################
######################################################################################################################
######################################################################################################################


###########################################################
# Load in PDX data
###########################################################

p_meta_t <- meta_t[meta_t$type == "PDX", ]
p_meta_k <- meta_k[meta_k$type == "PDX", ]

# load in arche zscores and normalize
zscore_t <- get_arche_scores("pdxs", "k50", p_meta_t)
normzs_t <- znorm(zscore_t)
zscore_k <- get_arche_scores("pdxs", "k50", p_meta_k)
normzs_k <- znorm(zscore_k)

# load in arche deviations and normalize
deviat_t <- get_arche_devs("PDXs", "50k", p_meta_t)
normdv_t <- znorm(deviat_t)
deviat_k <- get_arche_devs("PDXs", "50k", p_meta_k)
normdv_k <- znorm(deviat_k)

# pdx drug response data
xeva <- get_xeva("full")

###########################################################
# Removing low deviation samples
###########################################################

# get samples above threshold
zscore_t_sumdev <- get_devs(zscore_t, "zscore_pdx_t", 100, 150)
zscore_k_sumdev <- get_devs(zscore_k, "zscore_pdx_k", 100, 150)
deviat_t_sumdev <- get_devs(deviat_t, "deviat_pdx_t", 0.3, 0.5)
deviat_k_sumdev <- get_devs(deviat_k, "deviat_pdx_k", 0.3, 0.5)

# normalize
normzs_t_sumdev <- znorm(zscore_t_sumdev)
normzs_k_sumdev <- znorm(zscore_k_sumdev)
normdv_t_sumdev <- znorm(deviat_t_sumdev)
normdv_k_sumdev <- znorm(deviat_k_sumdev)

###########################################################
# Assign ARCHE scores to PDX xeva
###########################################################

# helper function
format_ARCHE <- function(xeva, scores) {
    scores <- t(scores) |> as.data.frame()
    df <- xeva[xeva$patient.id %in% rownames(scores),]
    df$ARCHE6 <- df$ARCHE5 <- df$ARCHE4 <- df$ARCHE3 <- df$ARCHE2 <- df$ARCHE1 <- NA
    for (i in 1:nrow(df)) {
        sample <- df$patient.id[i]
        df$ARCHE1[i] <- scores[rownames(scores) == sample,]$ARCHE1
        df$ARCHE2[i] <- scores[rownames(scores) == sample,]$ARCHE2
        df$ARCHE3[i] <- scores[rownames(scores) == sample,]$ARCHE3
        df$ARCHE4[i] <- scores[rownames(scores) == sample,]$ARCHE4
        df$ARCHE5[i] <- scores[rownames(scores) == sample,]$ARCHE5
        df$ARCHE6[i] <- scores[rownames(scores) == sample,]$ARCHE6
        
    }
    return(df)
}

# get ARCHE scores
zscore_t <- format_ARCHE(xeva, zscore_t)
normzs_t <- format_ARCHE(xeva, normzs_t)
zscore_k <- format_ARCHE(xeva, zscore_k)
normzs_k <- format_ARCHE(xeva, normzs_k)

deviat_t <- format_ARCHE(xeva, deviat_t)
normdv_t <- format_ARCHE(xeva, normdv_t)
deviat_k <- format_ARCHE(xeva, deviat_k)
normdv_k <- format_ARCHE(xeva, normdv_k)

zscore_t_sumdev <- format_ARCHE(xeva, zscore_t_sumdev)
normzs_t_sumdev <- format_ARCHE(xeva, normzs_t_sumdev)
zscore_k_sumdev <- format_ARCHE(xeva, zscore_k_sumdev)
normzs_k_sumdev <- format_ARCHE(xeva, normzs_k_sumdev)

deviat_t_sumdev <- format_ARCHE(xeva, deviat_t_sumdev)
normdv_t_sumdev <- format_ARCHE(xeva, normdv_t_sumdev)
deviat_k_sumdev <- format_ARCHE(xeva, deviat_k_sumdev)
normdv_k_sumdev <- format_ARCHE(xeva, normdv_k_sumdev)

###########################################################
# Compute PC of ARCHE-drug associations in PDXs
###########################################################

dir <- "data/results/figures/Misc/ARCHE_dr_test/"
x_zscore_t <- assess_ARCHE_PDX(zscore_t, "zscore_T", dir, plot = FALSE) |> suppressWarnings()
x_normzs_t <- assess_ARCHE_PDX(normzs_t, "normzscr_T", dir, plot = FALSE) |> suppressWarnings()
x_zscore_k <- assess_ARCHE_PDX(zscore_k, "zscore_K", dir, plot = FALSE) |> suppressWarnings()
x_normzs_k <- assess_ARCHE_PDX(normzs_k, "normzscr_K", dir, plot = FALSE) |> suppressWarnings()

x_deviat_t <- assess_ARCHE_PDX(deviat_t, "deviat_T", dir, plot = FALSE) |> suppressWarnings()
x_normdv_t <- assess_ARCHE_PDX(normdv_t, "normdev_T", dir, plot = FALSE) |> suppressWarnings()
x_deviat_k <- assess_ARCHE_PDX(deviat_k, "deviat_K", dir, plot = FALSE) |> suppressWarnings()
x_normdv_k <- assess_ARCHE_PDX(normdv_k, "normdev_K", dir, plot = FALSE) |> suppressWarnings()



x_zscore_t_sumdev <- assess_ARCHE_PDX(zscore_t_sumdev, "zscore_T", dir, plot = FALSE) |> suppressWarnings()
x_normzs_t_sumdev <- assess_ARCHE_PDX(normzs_t_sumdev, "normzscr_T", dir, plot = FALSE) |> suppressWarnings()
x_zscore_k_sumdev <- assess_ARCHE_PDX(zscore_k_sumdev, "zscore_K", dir, plot = FALSE) |> suppressWarnings()
x_normzs_k_sumdev <- assess_ARCHE_PDX(normzs_k_sumdev, "normzscr_K", dir, plot = FALSE) |> suppressWarnings()

x_deviat_t_sumdev <- assess_ARCHE_PDX(deviat_t_sumdev, "deviat_T", dir, plot = FALSE) |> suppressWarnings()
x_normdv_t_sumdev <- assess_ARCHE_PDX(normdv_t_sumdev, "normdev_T", dir, plot = FALSE) |> suppressWarnings()
x_deviat_k_sumdev <- assess_ARCHE_PDX(deviat_k_sumdev, "deviat_K", dir, plot = FALSE) |> suppressWarnings()
x_normdv_k_sumdev <- assess_ARCHE_PDX(normdv_k_sumdev, "normdev_K", dir, plot = FALSE) |> suppressWarnings()

###########################################################
# Compile PDX results
###########################################################

# all samples
pdx_all_toPlot <- rbind(
    x_zscore_t, x_normzs_t,
    x_zscore_k, x_normzs_k,
    x_deviat_t, x_normdv_t,
    x_deviat_k, x_normdv_k
)
colnames(pdx_all_toPlot)[colnames(pdx_all_toPlot) == "ARCHE_label"] <- "Label"
colnames(pdx_all_toPlot)[colnames(pdx_all_toPlot) == "pair"] <- "pairs"

# subsetted samples
pdx_sumdev_toPlot <- rbind(
    x_zscore_t_sumdev, x_normzs_t_sumdev,
    x_zscore_k_sumdev, x_normzs_k_sumdev,
    x_deviat_t_sumdev, x_normdv_t_sumdev,
    x_deviat_k_sumdev, x_normdv_k_sumdev
)
colnames(pdx_sumdev_toPlot)[colnames(pdx_sumdev_toPlot) == "ARCHE_label"] <- "Label"
colnames(pdx_sumdev_toPlot)[colnames(pdx_sumdev_toPlot) == "pair"] <- "pairs"


# last save: April 22, 2026
save(cell_toPlot, pdx_all_toPlot, pdx_sumdev_toPlot, file = "data/results/data/Misc/ARCHE_drug_response_testing.RData")

###########################################################
# Plot PDX results
###########################################################

# helper function to plot PDX ARCHE drug response associations
plot_PDX_ARCHEs <- function(compile, label) {

    sig <- compile[which(abs(compile$PC.BAR_median) > 0.4 & compile$pval.BAR_median < 0.1),]
    sig_pairs <- sig$pair
    toPlot <- compile[compile$pair %in% sig_pairs,]
    toPlot$sig <- ifelse(toPlot$pval.BAR_median < 0.1, 'pval < 0.1', 'pval >= 0.1')

    toPlot$Label <- factor(toPlot$Label, levels = c(
        "zscore_T", "zscore_K", "deviat_T", "deviat_K",
        "normzscr_T", "normzscr_K", "normdev_T", "normdev_K")
    )

    for (arche in paste0("ARCHE", 1:6)) {

        subset <- toPlot[toPlot$ARCHE == arche,]
        subset$score <- ifelse(
            sub("_.*", "", subset$Label) %in% c("zscore", "normzscr"),
            ifelse(sub(".*_", "", subset$Label) == "T", "ZScore_T", "ZScore_K"),
            ifelse(sub(".*_", "", subset$Label) == "T", "Dev_T", "Dev_K"))
        

        p1 <- ggplot(subset, aes(x = drug, y = score, fill = N)) +
            geom_tile() +
            geom_text(aes(label = N), size = 2.5) +
            scale_fill_gradient(high = "#B6B8D6", low = "#BBDBD1") +
            theme_void() +
            theme(
                axis.text.y = element_text(size = 8, hjust = 1),
                legend.position = "none") +
            ggtitle(paste0(arche, "_", label))

        p2 <- ggplot(subset, aes(x = drug, y = Label, fill = PC.BAR_median, size = -log(pval.BAR_median), shape = sig)) +
            geom_point() +
            geom_text(data = subset(subset, sig == 'pval < 0.1'), aes(label = round(PC.BAR_median, 2)), color = "black", size = 2.5) +
            scale_shape_manual(values = c(21, 24)) +
            scale_size(range = c(2, 12)) +
            scale_fill_gradient2(
                low = "#BC4749",
                high = "#689CB0",
                mid = "#C2BBC9",
                limits = c(-1, 1)
            ) +
            theme_bw() +
            theme(
                legend.key.size = unit(0.3, 'cm'),
                axis.text.x = element_text(size=6, angle=25, hjust=1, vjust=1, margin = margin(t = 3))
            ) 
        
        p <- p1 / p2 + plot_layout(heights = c(1, 6))
        filename <- paste0("data/results/figures/Misc/ARCHE_dr_test/PDXs_", arche, "_", label, ".png")
        ggsave(filename, p, w = 10, h = 5)
    }

}

plot_PDX_ARCHEs(pdx_all_toPlot, "all_samples")
plot_PDX_ARCHEs(pdx_sumdev_toPlot, "sumdev_samples")
