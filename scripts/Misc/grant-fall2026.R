# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(patchwork)
    library(ggnewscale)
    library(ggpattern)
})

source("utils/get_data.R")
source("utils/palettes.R")
source("utils/mappings.R")
source("utils/plots/ARCHE_scores_heatmap.R")
source("utils/ccl_benchmarks.R")

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
zscore_cells <- get_arche_scores("data/rawdata/all_scoring/cell_tcga.Zscore.txt", c_meta)
normzs_cells <- znorm(zscore_cells)

# sumdevs
zscore_cells_sumdev <- get_arche_sumdevs(zscore_cells, "zscore_cells", plot = TRUE)
normzs_cells_sumdev <- znorm(zscore_cells_sumdev)

# load in RNA
ubr1 <- get_pset_rna("UBR1")
ubr2 <- get_pset_rna("UBR2")
gray <- get_pset_rna("GRAY")
gcsi <- get_pset_rna("gCSI")
gdsc <- get_pset_rna("GDSC2")
ccle <- get_pset_rna("CCLE")

# load in mutation data
gdsc_mut <- get_pset_mut("GDSC2")
ccle_mut <- get_pset_mut("CCLE")

# load in PAM50 subtyping
load("data/results/data/3-DataExploration/ccls_subtyping_scores.RData")

# get drug sensitivity data
load("data/procdata/CCLs/sensitivity_data.RData")

# load in PCC
load("data/results/data/4-DrugResponse/CCLs/ARCHE_CCLs_associations.RData")


###########################################################
# Compile PAM50 data
###########################################################

pam50_scores <- rbind(
    format_pam50(ubr1_pam50, "UBR1"), format_pam50(ubr2_pam50, "UBR2"), format_pam50(gray_pam50, "GRAY"),
    format_pam50(gcsi_pam50, "gCSI"), format_pam50(gdsc_pam50, "GDSC2"), format_pam50(ccle_pam50, "CCLE"),
    format_pam50(ccle_pam50, "CTRP")
)


###########################################################
# Helper functions
###########################################################

plot_ccl_benchmark <- function(arche, ft_comp, drug, ft_mat, ft_lab, arche_scores, pset_arche, pset_comp) {

    comp_lab <- ifelse(ft_lab == "RNA", "ARCHE", "PAM50")

    # get individual associations
    arche_sen <- get_all_drug_sen(paste0(arche, "_", drug), arche_scores, "ARCHE", c_meta)
    comp_sen <- get_all_drug_sen(paste0(ft_comp, "_", drug), ft_mat, comp_lab, c_meta)

    # plot ARCHE and comparison feature
    p1 <- plot_scatter(paste0(arche, "_", drug), arche_sen, pset_arche, "ARCHE") + theme(legend.position = "none")
    p2 <- plot_scatter(paste0(ft_comp, "_", drug), comp_sen, pset_comp, ft_lab) + theme(axis.title.y = element_blank())
    p <- p1 + p2
    filename <- paste0("data/results/figures/Misc/cihr-fall2026/", arche, "_", drug, ".png")
    ggsave(filename, p, width = 7, height = 2.5)

}

plot_scatter <- function(pair, all_drug_sen, pset, label) {

    all_drug_sen <- all_drug_sen[all_drug_sen$PSet == pset,]

    corr <- cor.test(all_drug_sen$AAC, all_drug_sen$Score, method = "pearson")
    pcc <- round(corr$estimate, 4)
    pval <- round(corr$p.value, 5)

    feature <- gsub("_.*", "", pair)
    drug <- gsub(paste0(feature, "_"), "", pair)
    x_label <- ifelse(label == "RNA", "Expression", "Score")

    all_drug_sen$PSet <- factor(all_drug_sen$PSet, levels = names(PSet_pal))
    p <- ggplot(all_drug_sen, aes(x = Score, y = AAC, fill = Subtype)) + 
        geom_point(size = 3, shape = 21) + 
        geom_smooth(method = "lm", se=TRUE, color = "black", aes(group = 1), show.legend = FALSE) + 
        scale_fill_manual(values = subtype_pal) +
        guides(fill = guide_legend(override.aes = list(size = 4))) +
        theme_bw() + 
        theme(
            strip.background = element_rect(fill = "white"),
            legend.key.size = unit(0.5, 'cm'),
            axis.title.y = element_text(size = 11, margin = margin(r = 10))
        ) +
        labs(
            x = paste(feature, x_label), 
            y = paste(drug, "Response (AAC)"),
            title = paste("PCC:", pcc, ", pval:", pval)
        )
    return(p)
}

plot_mut_scatter <- function(arche, drug, genes, arche_scores, pset, mut_df) {

    arche_sen <- get_all_drug_sen(paste0(arche, "_", drug), arche_scores, "ARCHE", c_meta)

    mut_df <- as.data.frame(t(mut_df[names(genes),]))
    subset <- arche_sen[arche_sen$PSet == pset & arche_sen$Sample %in% rownames(mut_df),]
    subset$Mut <- rowSums(mut_df)[match(subset$Sample, names(rowSums(mut_df)))]
    subset <- subset[!is.na(subset$Mut),]
    subset$Mut <- ifelse(subset$Mut > 0, 1, 0)
    for (gene in names(genes)) {
        subset[[gene]] <- mut_df[[gene]][match(subset$Sample, rownames(mut_df))]
    }

    subset$Mut <- factor(subset$Mut, levels = c(1, 0), labels = c("Mut", "WT"))
    subset <- subset[order(subset$AAC, decreasing = TRUE),]
    subset$Sample <- factor(subset$Sample, levels = unique(subset$Sample))
    subset$dummy <- c("BRAF", "NF1", rep("WT", nrow(subset)-2))

    corr <- cor.test(subset$AAC, subset$Score, method = "pearson")
    estimate <- round(corr$estimate, 4)
    pval <- round(corr$p.value, 5)
    subset$Subtype <- factor(subset$Subtype, levels = names(subtype_pal))

    p <- ggplot(subset, aes(x = Score, y = AAC, fill = Subtype)) +
        geom_smooth(method = "lm", se = TRUE, color = "black", aes(group = 1), show.legend = FALSE) +
        geom_point(size = 3.5, shape = 21) +
        scale_fill_manual(values = subtype_pal) +
        new_scale_color() +
        geom_point(data = subset[subset$BRAF == 1,], aes(x = Score, y = AAC), size = 5.5, shape = 0) +
        geom_point(data = subset[subset$NF1 == 1,], aes(x = Score, y = AAC), size = 5.5, shape = 5) +
        # for the legend
        geom_point(data = subset, aes(x = Score, y = AAC, shape = dummy), size = 0, alpha = 0) +
        scale_shape_manual("Mutation", values = c("BRAF" = 0, "NF1" = 5, "WT" = 1)) +
        guides(
            fill = guide_legend(override.aes = list(size = 4, shape = 21)),
            shape = guide_legend(override.aes = list(size = 4, alpha = 1))
        ) +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.title.y = element_text(size = 11, margin = margin(r = 10)),
            axis.title.x = element_text(size = 11)
        ) +
        labs(y = paste(drug, "Response (AAC)"), x = paste(arche, "Score"), title = paste0("PCC: ", estimate, " | pval: ", pval))
    filename <- paste0("data/results/figures/Misc/cihr-fall2026/", arche, "_", drug, ".png")
    ggsave(filename, p, width = 3.2, height = 2.5)
}


###########################################################
# Plot associations
###########################################################

plot_ccl_benchmark("ARCHE5", "Basal", "Paclitaxel", pam50_scores, "PAM50", normzs_cells_sumdev, "GRAY", "GRAY")
plot_ccl_benchmark("ARCHE4", "ENSG00000198900", "Topotecan", ccle, "RNA", normzs_cells_sumdev, "GDSC2", "CCLE")

genes <- c("NF1" = "ENSG00000196712", "BRAF" = "ENSG00000157764")

plot_mut_scatter("ARCHE4", "Selumetinib", genes, normzs_cells_sumdev, "GDSC2", gdsc_mut)
plot_mut_scatter("ARCHE4", "Trametinib", genes, normzs_cells_sumdev, "GDSC2", gdsc_mut)
plot_mut_scatter("ARCHE4", "PD-0325901", genes, normzs_cells_sumdev, "GDSC2", gdsc_mut)