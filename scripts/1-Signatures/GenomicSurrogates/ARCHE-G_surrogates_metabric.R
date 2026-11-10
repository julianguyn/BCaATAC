suppressPackageStartupMessages({
    library(data.table)
    library(ComplexHeatmap)
    library(circlize)
    library(survival)
    library(ggsurvfit)
    library(survminer)
    library(ggplot2)
    library(dplyr)
    library(patchwork)
    library(genefu)
    library(org.Hs.eg.db)
})

source("utils/get_data.R")
source("utils/palettes.R")
source("utils/gene_surrogate/helpers.R")

set.seed(123)

dir <- "data/results/data/1-Signatures/geneSurrogates"
fig_dir <- "data/results/figures/1-Signatures/metabric"
dir.create(paste0(fig_dir, "/ARCHEG_survivalplots"), recursive = TRUE, showWarnings = FALSE)

###########################################################
# Plot model accuracy from training (TCGA CV)
###########################################################

# per-fold performance across 10 repeats of 5-fold CV
folds <- data.frame(matrix(nrow=0, ncol=0))
for (file in list.files(dir, pattern = "en_folds.csv")) {
    df <- read.csv(paste0(dir, "/", file))
    df$ARCHE <- toupper(sub("_.*", "", file))
    folds <- rbind(folds, df)
}

p <- ggplot(folds, aes(x = ARCHE, y = spearman, fill = ARCHE)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
    scale_fill_manual(values = ARCHE_pal) +
    labs(x = NULL, y = "Spearman (held-out fold)") +
    theme_classic(base_size = 13) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("data/results/figures/1-Signatures/ARCHE-G_training_accuracy.png", p, width = 4.5, height = 3.5)

# observed vs predicted, averaging each sample's prediction across repeats
cv_preds <- data.frame(matrix(nrow=0, ncol=0))
for (file in list.files(dir, pattern = "en_predictions.csv")) {
    df <- read.csv(paste0(dir, "/", file))
    df <- df %>% group_by(sample) %>% summarise(y_true = dplyr::first(y_true), y_pred = mean(y_pred))
    df$ARCHE <- toupper(sub("_.*", "", file))
    cv_preds <- rbind(cv_preds, df)
}

labels <- cv_preds %>%
    group_by(ARCHE) %>%
    summarise(
        rho = cor(y_true, y_pred, method = "spearman"),
        R2 = 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)
    ) %>%
    mutate(label = sprintf("rho = %.2f\nR2 = %.2f", rho, R2))

p <- ggplot(cv_preds, aes(x = y_true, y = y_pred, color = ARCHE)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_text(data = labels, aes(x = -Inf, y = Inf, label = label), hjust = -0.1, vjust = 1.2, size = 3, color = "black") +
    facet_wrap(~ARCHE, nrow = 2) +
    scale_color_manual(values = ARCHE_pal) +
    labs(x = "Observed ARCHE score", y = "CV-predicted ARCHE score") +
    theme_classic(base_size = 13) +
    theme(legend.position = "none")
ggsave("data/results/figures/1-Signatures/ARCHE-G_training_obs_vs_pred.png", p, width = 7, height = 5)

###########################################################
# Load in METABRIC data
###########################################################

# clinical data
pheno <- read.table("data/rawdata/metabric/data_clinical_patient.txt", header = TRUE, sep = "\t", comment.char = "#", quote = "", na.strings = c("", "NA"))
rownames(pheno) <- gsub("-", "\\.", pheno$PATIENT_ID)
pheno$overall_survival <- pheno$OS_MONTHS * 30.4375 # convert to days to match tcga
pheno$status <- as.numeric(sub(":.*", "", pheno$OS_STATUS))

# rna
rna <- fread("data/procdata/metabric/METABRIC_BRCA_Illumina_HT12v3_microarray_log2.cct", data.table = FALSE)
rownames(rna) <- rna$attrib_name
rna <- t(rna[,-1])

# fill few missing values with the gene mean
na_idx <- which(is.na(rna), arr.ind = TRUE)
rna[na_idx] <- colMeans(rna, na.rm = TRUE)[na_idx[,2]]

# load in ARCHE-G predictions
preds <- read.csv(paste0(dir, "/metabric/METABRIC_BRCA_Illumina_HT12v3_microarray_log2_ARCHE_predictions.csv"), row.names = 1)

###########################################################
# PAM50 subtyping
###########################################################

# load in genefu subtyping model
data(pam50.robust)

entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = colnames(rna),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
)

annot <- data.frame(
    probe = colnames(rna),
    EntrezGene.ID = as.integer(entrez_ids),
    row.names = colnames(rna),
    stringsAsFactors = FALSE
)

pam50 <- molecular.subtyping(
    sbt.model = "pam50",
    data = rna,
    annot = annot,
    do.mapping = TRUE
)
saveRDS(pam50, file = "data/procdata/metabric/pam50_subtyping.rds")

pam50_scores <- as.data.frame(pam50$subtype.proba)
pam50_scores$assigned <- as.character(pam50$subtype)

# compare with cbioportal pam50 + claudin-low calls
table(genefu = pam50_scores$assigned, cbioportal = pheno$CLAUDIN_SUBTYPE[match(rownames(pam50_scores), rownames(pheno))])

###########################################################
# Correlate ARCHE-G and PAM50 scores
###########################################################

scores <- pam50_scores[match(rownames(preds), rownames(pam50_scores)),]

compiled <- data.frame(matrix(nrow=0, ncol=0))
for (arche in colnames(preds)) {
    tobind <- data.frame(
        Basal = cor.test(preds[[arche]], scores$Basal)$estimate,
        Her2 = cor.test(preds[[arche]], scores$Her2)$estimate,
        LumA = cor.test(preds[[arche]], scores$LumA)$estimate,
        LumB = cor.test(preds[[arche]], scores$LumB)$estimate,
        Normal = cor.test(preds[[arche]], scores$Normal)$estimate
    )
    rownames(tobind) <- arche
    compiled <- rbind(compiled, tobind)
}

###########################################################
# Plot compiled correlations
###########################################################

toPlot <- t(compiled)

cols <- colorRampPalette(c("#39066B", "#B18ED6", "#EBEBEB", "#75BFB2", "#008A8C"))(9)
col_fun <- colorRamp2(seq(-1,1,length.out = 9),cols)

col_ha <- HeatmapAnnotation(
    'ARCHE' = paste0("ARCHE", 1:6),
    col = list('ARCHE' = ARCHE_pal),
    simple_anno_size = unit(3, "mm"),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

row_ha <- rowAnnotation(
    Subtype = names(subtype_pal)[1:5],
    col = list(Subtype = subtype_pal),
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size = unit(3, "mm"),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

ht <- Heatmap(
    toPlot,
    column_title_gp = gpar(fontsize = 9),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "Correlation",
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 8),
    rect_gp = gpar(col = "white", lwd = 0.5),
    bottom_annotation = col_ha,
    left_annotation = row_ha,
)
filename <- paste0(fig_dir, "/ARCHE-G_subtype_heatmap.png")
png(filename, width = 3, height = 2, res = 600, units = "in")
ht
dev.off()

###########################################################
# Assign ARCHE-G
###########################################################

archeG <- data.frame(
    Sample = rownames(preds),
    ARCHEG = colnames(preds)[max.col(preds, ties.method = "first")],
    Score = apply(preds, 1, max)
)
archeG$PAM50 <- pam50_scores$assigned[match(archeG$Sample, rownames(pam50_scores))]
saveRDS(archeG, file = "data/procdata/metabric/ARCHEG_metabric.rds")

###########################################################
# Plot ARCHE-G PAM50 subtype counts
###########################################################

toPlot <- t(as.matrix(table(archeG$ARCHEG, archeG$PAM50)))

cols <- colorRampPalette(c("#EBEBEB", "#75BFB2", "#008A8C"))(9)
col_fun <- colorRamp2(seq(0,max(toPlot),length.out = 9),cols)

ht <- Heatmap(
    toPlot,
    column_title_gp = gpar(fontsize = 9),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "No. Samples",
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 8),
    rect_gp = gpar(col = "white", lwd = 0.5),
    bottom_annotation = col_ha,
    left_annotation = row_ha,
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (toPlot[i, j] > 0) {
            grid.text(toPlot[i, j], x, y, gp = gpar(fontsize = 7, col = "black"))
        }
    }
)
filename <- paste0(fig_dir, "/ARCHE-G_subtype_count_heatmap.png")
png(filename, width = 3, height = 2, res = 600, units = "in")
ht
dev.off()

###########################################################
# Plot survival plots
###########################################################

# plot both 5 and 10 years
HR_5 <- plot_survival(archeG, pheno, 5, outdir = paste0(fig_dir, "/ARCHEG_survivalplots"), highlight = FALSE)
HR_10 <- plot_survival(archeG, pheno, 10, outdir = paste0(fig_dir, "/ARCHEG_survivalplots"), highlight = FALSE)
