# score mammary lineage, pluripotency, WNT and neural crest signatures (ssGSEA) across ARCHEs

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(GSVA)
    library(msigdbr)
    library(ComplexHeatmap)
    library(circlize)
})

source("utils/get_data.R")
source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# ARCHE labels for all tumours
arche_meta <- read.csv("data/rawdata/tcga/TCGA_sourcefiles.csv")

# load in tumour gene counts matrix (log2 RSEM)
t_counts <- get_tcga_rna()

###########################################################
# Get MSigDB signatures
###########################################################

# signatures grouped by the question each figure answers
# (WP_WNT_SIGNALING_PATHWAY is named WP_WNT_SIGNALING in current MSigDB)
signature_groups <- list(
    mammary_lineage = c(
        "PECE_MAMMARY_STEM_CELL_UP", "PECE_MAMMARY_STEM_CELL_DN",
        "LIM_MAMMARY_STEM_CELL_UP", "LIM_MAMMARY_STEM_CELL_DN",
        "LIM_MAMMARY_LUMINAL_PROGENITOR_UP",
        "LIM_MAMMARY_LUMINAL_MATURE_UP", "LIM_MAMMARY_LUMINAL_MATURE_DN"
    ),
    pluripotency = c(
        "BENPORATH_ES_1", "BENPORATH_NOS_TARGETS", "BENPORATH_PRC2_TARGETS",
        "WONG_EMBRYONIC_STEM_CELL_CORE", "MUELLER_PLURINET"
    ),
    wnt = c(
        "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY", "GOBP_NON_CANONICAL_WNT_SIGNALING_PATHWAY", "WP_WNT_SIGNALING"
    ),
    neural_crest = c(
        "GOBP_NEURAL_CREST_CELL_DIFFERENTIATION", "GOBP_NEURAL_CREST_CELL_FATE_COMMITMENT",
        "WP_NEURAL_CREST_DIFFERENTIATION", "WP_NEURAL_CREST_CELL_MIGRATION_DURING_DEVELOPMENT",
        "LEE_NEURAL_CREST_STEM_CELL_UP", "LEE_NEURAL_CREST_STEM_CELL_DN"
    ),
    great_rna_checks = c(
        "GOBP_METANEPHROS_DEVELOPMENT", "GOBP_CILIUM_MOVEMENT",
        "TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_DUCTAL_NORMAL_UP"
    )
)
signatures <- unlist(signature_groups, use.names = FALSE)

msig <- rbind(
    msigdbr(species = "Homo sapiens", collection = "C2"),
    msigdbr(species = "Homo sapiens", collection = "C5")
)
msig <- msig[msig$gs_name %in% signatures,]
gene_sets <- lapply(split(msig$gene_symbol, msig$gs_name), unique) # remove duplicate genes per set

###########################################################
# Single sample gsea
###########################################################

# ssgsea on tumour counts (signatures x samples)
t_es <- gsva(ssgseaParam(t_counts, gene_sets), verbose = FALSE)

# z-score each signature across tumours so signatures are on the same scale
t_es <- t(scale(t(t_es)))

###########################################################
# Plot signatures by ARCHE
###########################################################

scores <- as.data.frame(t(t_es))
scores$ARCHE <- arche_meta$ARCHE[match(rownames(scores), arche_meta$Sample.Name)]
scores <- scores[!is.na(scores$ARCHE),]

# mean score of tumours in each ARCHE (signatures x ARCHEs)
arche_means <- aggregate(. ~ ARCHE, data = scores, FUN = mean)
rownames(arche_means) <- arche_means$ARCHE
arche_means <- t(as.matrix(arche_means[paste0("ARCHE", 1:6), -1]))

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#B2182B", "#F7F7F7", "#2166AC"))

ha <- HeatmapAnnotation(
    'ARCHE' = colnames(arche_means),
    col = list('ARCHE' = ARCHE_pal),
    simple_anno_size = unit(3, "mm"),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

# one heatmap per group of signatures
# names(signature_groups)
plot_heatmap <- function(group, w, h) {
    toPlot <- arche_means[signature_groups[[group]], , drop = FALSE]

    ht <- Heatmap(
        toPlot,
        name = "Mean\nZ-score",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8),
        row_names_max_width = max_text_width(rownames(toPlot), gp = gpar(fontsize = 8)), # fit long signature names
        column_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 0.5),
        bottom_annotation = ha
    )

    filename <- paste0("data/results/figures/2-MolecularSigAnalysis/molecularsig/stemness_", group, "_heatmap.png")
    png(filename, width = w, height = h, res = 600, units = "in")
    draw(ht)
    dev.off()
}

plot_heatmap("mammary_lineage", 4, 2)
plot_heatmap("pluripotency", 3.75, 1.5)
plot_heatmap("wnt", 4.5, 1.25)
plot_heatmap("neural_crest", 5, 1.6)
plot_heatmap("great_rna_checks", 5.5, 1.25)
