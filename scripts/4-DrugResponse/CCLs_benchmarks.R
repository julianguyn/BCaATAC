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
source("utils/compute_drug_response.R")
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

# sumdev zscores
zscore_cells <- get_arche_scores(paste0("data/rawdata/all_scoring/cell_tcga.Zscore.txt"), c_meta)
zscore_cells_sumdev <- get_arche_sumdevs(zscore_cells)

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

###########################################################
# Compile PAM50 data
###########################################################

pam50_scores <- rbind(
    format_pam50(ubr1_pam50, "UBR1"), format_pam50(ubr2_pam50, "UBR2"), format_pam50(gray_pam50, "GRAY"),
    format_pam50(gcsi_pam50, "gCSI"), format_pam50(gdsc_pam50, "GDSC2"), format_pam50(ccle_pam50, "CCLE"),
    format_pam50(ccle_pam50, "CTRP")
)

###########################################################
# Plot associations
###########################################################

tozasertib_genes <- c(
    "ENSG00000087586" = "AURKA",
    "ENSG00000178999" = "AURKB",
    "ENSG00000105146" = "AURKC"
)

paclitaxel_genes <- c(
    "ENSG00000258947" = "TUBB3",
    "ENSG00000117632" = "STMN1",
    "ENSG00000047849" = "MAP4",
    "ENSG00000085563" = "ABCB1"
)

trastuzumab_genes <- c(
    "ENSG00000141736" = "ERBB2"
)

erlotinib_genes <- c(
    "ENSG00000146648" = "EGFR",
    "ENSG00000123374" = "CDK2"
)

topotecan_genes <- c(
    "ENSG00000172716" = "SLFN11"
)

PD_0325901_genes <- c(
    "ENSG00000157764" = "BRAF",
    "ENSG00000213281" = "NRAS",
    "ENSG00000133703" = "KRAS",
    "ENSG00000196712" = "NF1",
    "ENSG00000169032" = "MAP2K1"
)

# ARCHE2
a2_tozasertib <- plot_rna_associations("ARCHE2", "Basal", "Tozasertib", zscore_cells_sumdev, tozasertib_genes)
a2_topotecan <- plot_rna_associations("ARCHE2", "Basal", "Topotecan", zscore_cells_sumdev, topotecan_genes)

# ARCHE3
a3_trastuzumab <- plot_rna_associations("ARCHE3", "Her2", "Trastuzumab", zscore_cells_sumdev, trastuzumab_genes)

# ARCHE4
a4_erlotinib <- plot_rna_associations("ARCHE4", "Basal", "Erlotinib", zscore_cells_sumdev, erlotinib_genes)
a4_topotecan <- plot_rna_associations("ARCHE4", "Basal", "Topotecan", zscore_cells_sumdev, topotecan_genes)
a4_topotecan <- plot_rna_associations("ARCHE4", "Basal", "Topotecan", zscore_cells_sumdev, topotecan_genes)
a4_PD_0325901 <- plot_rna_associations("ARCHE4", "LumA", "PD-0325901", zscore_cells_sumdev, PD_0325901_genes)
a4_trametinib <- plot_rna_associations("ARCHE4", "Basal", "Trametinib", zscore_cells_sumdev, PD_0325901_genes)


# ARCHE5
a5_paclitaxel <- plot_rna_associations("ARCHE5", "Basal", "Paclitaxel", zscore_cells_sumdev, paclitaxel_genes)
a5_erlotinib <- plot_rna_associations("ARCHE5", "Basal", "Erlotinib", zscore_cells_sumdev, erlotinib_genes)


###########################################################
# Mutations
###########################################################

# these won't work anymore
#plot_mut_associations("ARCHE3", "Alpelisib", "PIK3CA", zscore_cells_sumdev)
#plot_mut_associations("ARCHE4", "Erlotinib", "EGFR", zscore_cells_sumdev)

mutations = c("NF1", "BRAF")

plot_mut_associations("ARCHE4", "Selumetinib", mutations, zscore_cells_sumdev, "GDSC2", "CCLE", "CTRP")
plot_mut_associations("ARCHE4", "PD-0325901", mutations, zscore_cells_sumdev, "gCSI", "GDSC2", "CCLE")
plot_mut_associations("ARCHE4", "Trametinib", mutations, zscore_cells_sumdev, "GRAY", "GDSC2", "CTRP")

#' Helper function
plot_mut <- function(arche, drug, mut, pset, arche_sen) {

    mut_df <- switch(
        pset,
        GDSC2 = gdsc_mut,
        CCLE = ccle_mut,
        CTRP = ccle_mut,
        GRAY = NA,
        gCSI = NA
    )

    # ---------- With mutation data
    if (class(mut_df) == "data.frame") {

        mut_df <- as.data.frame(t(mut_df[mutations,]))
        subset <- arche_sen[arche_sen$PSet == pset & arche_sen$Sample %in% rownames(mut_df),]
        subset$Mut <- rowSums(mut_df)[match(subset$Sample, names(rowSums(mut_df)))]
        subset <- subset[!is.na(subset$Mut),]
        subset$Mut <- ifelse(subset$Mut > 0, 1, 0)
        for (mut in mutations) {
            subset[[mut]] <- mut_df[[mut]][match(subset$Sample, rownames(mut_df))]
        }

        subset$Mut <- factor(subset$Mut, levels = c(1, 0), labels = c("Mut", "WT"))
        subset <- subset[order(subset$AAC, decreasing = TRUE),]
        subset$Sample <- factor(subset$Sample, levels = unique(subset$Sample))
        subset$dummy <- c("BRAF", "NF1", rep("WT", nrow(subset)-2))

        p1 <- ggplot(subset, aes(x = Sample, y = AAC, pattern = Mut, fill = Subtype)) +
            geom_col_pattern(
                alpha = 0.85, pattern_fill = "black",
                pattern_density = 0.5, pattern_spacing = 0.03,
                pattern_color = NA, color = "black"
            ) +
            scale_pattern_manual(values = c(WT = "none", Mut = "stripe")) +
            scale_fill_manual(values = subtype_pal) +
            theme_void() +
            theme(
                axis.text.x = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 11),
                axis.title.y = element_blank(), axis.title.x = element_text(size = 9),
                legend.position = "none"
            ) + 
            labs(title = pset, x = "Samples ranked by AAC")

        p2 <- ggplot(subset, aes(x = Score, y = AAC, fill = Subtype)) +
            geom_smooth(method = "lm", se = TRUE, color = "black", aes(group = 1), show.legend = FALSE) +
            geom_point(size = 2.5, shape = 21) +
            scale_fill_manual(values = subtype_pal) +
            new_scale_color() +
            geom_point(data = subset[subset$BRAF == 1,], aes(x = Score, y = AAC), size = 4, shape = 0) +
            geom_point(data = subset[subset$NF1 == 1,], aes(x = Score, y = AAC), size = 4, shape = 5) +
            # for the legend
            geom_point(data = subset, aes(x = Score, y = AAC, shape = dummy), size = 0, alpha = 0) +
            scale_shape_manual("Mutation", values = c("BRAF" = 0, "NF1" = 5, "WT" = 1)) +
            guides(
                fill = guide_legend(override.aes = list(size = 4, shape = 21)),
                shape = guide_legend(override.aes = list(size = 4, alpha = 1))
            ) +
            theme_bw() +
            theme(
                legend.key = element_blank(),
                legend.key.size = unit(0.5, 'cm'),
                axis.title.y = element_text(size = 9, margin = margin(r = 10)),
                axis.title.x = element_text(size = 9)
            ) +
            ylim(0, 0.5) +
            labs(y = paste(drug, "Response (AAC)"), x = paste(arche, "Score"))


    } else { # ---------- Without mutation data

        subset <- arche_sen[arche_sen$PSet == pset,]
        subset$Mut <- c("Mut", rep("WT", nrow(subset)-1))

        p1 <- ggplot(subset, aes(x = Sample, y = AAC, pattern = Mut, fill = Subtype)) +
            geom_col_pattern(
                alpha = 0.85, pattern_fill = "black",
                pattern_density = 0.5, pattern_spacing = 0.03,
                pattern_color = NA, color = "black"
            ) +
            scale_pattern_manual(values = c(WT = "none", Mut = "stripe")) +
            scale_fill_manual(values = subtype_pal) +
            theme_void() +
            theme(
                axis.text.x = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 11),
                axis.title.y = element_blank(), axis.title.x = element_text(size = 9),
                legend.position = "none"
            ) + 
            labs(title = pset, x = "Samples ranked by AAC")
        p2 <- ggplot(subset, aes(x = Score, y = AAC, fill = Subtype)) +
            geom_smooth(method = "lm", se = TRUE, color = "black", aes(group = 1), show.legend = FALSE) +
            geom_point(size = 2.5, shape = 21) +
            scale_fill_manual(values = subtype_pal) +
            theme_bw() +
            theme(
                legend.key = element_blank(),
                legend.key.size = unit(0.5, 'cm'),
                axis.title.y = element_text(size = 9, margin = margin(r = 10)),
                axis.title.x = element_text(size = 9)
            ) +
            ylim(0, 0.5) +
            labs(y = paste(drug, "Response (AAC)"), x = paste(arche, "Score"))
    }

    p <- p1 / p2 + plot_layout(height = c(1, 2))
    return(p)
}


#' Plot mutation associations
plot_mut_associations <- function(arche, drug, mut, arche_scores, pset1, pset2, pset3) {
    arche_sen <- get_all_drug_sen(paste0(arche, "_", drug), arche_scores, "ARCHE", c_meta)

    p1 <- plot_mut(arche, drug, mut, pset1, arche_sen) + theme(legend.position = "none")
    p2 <- plot_mut(arche, drug, mut, pset2, arche_sen) + theme(legend.position = "none") + theme(axis.title.y = element_blank())
    p3 <- plot_mut(arche, drug, mut, pset3, arche_sen) + theme(axis.title.y = element_blank())

    p <- wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3) + 
        plot_layout(width = c(1.05,1,1.3))
    filename <- paste0("data/results/figures/4-DrugResponse/benchmarks/", arche, "_", drug, "/mut_compiled.png")
    ggsave(filename, p, width = 10, height = 3.75)
}
