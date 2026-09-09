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
})

source("utils/get_data.R")
source("utils/palettes.R")
source("utils/gene_surrogate/helpers.R")

set.seed(123)

dir <- "data/results/data/1-Signatures/geneSurrogates"

###########################################################
# Load in data
###########################################################

# read in tcga clinical data
pheno <- read.table("data/rawdata/tcga/Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", header = TRUE, row.names = 1)
pheno <- t(pheno) |> as.data.frame()
pheno$overall_survival <- as.numeric(pheno$overall_survival)
pheno$status <- as.numeric(pheno$status)

# load in pam50 subtypes
pam50 <- readRDS("data/procdata/TCGA/pam50_subtyping_full_cohort.rds")
pam50_scores <- as.data.frame(pam50$subtype.proba)
pam50_scores$assigned <- as.character(pam50$subtype)

###########################################################
# Load in CV predictions
###########################################################

compiled <- data.frame(matrix(nrow=0, ncol=0))

R2 <- c()

for (file in list.files(dir, pattern = "en_predictions.csv")) {
    arche <- toupper(sub("_.*", "", file))
    df <- read.csv(paste0(dir, "/", file))
    df$ARCHE <- arche
    df$error <- df$y_true - df$y_pred

    errors <- c()
    y_true <- c()

    for (sample in unique(df$sample)) {
        subset_df <- df[df$sample == sample,]
        errors <- c(errors, subset_df$error[abs(subset_df$error) == min(abs(subset_df$error))])
        y_true <- c(y_true, subset_df$y_true[subset_df$error == min(subset_df$error)])
    }
    R2 <- c(R2, 1 - sum(errors^2) / sum((y_true - mean(y_true))^2))
}

###########################################################
# Correlate ARCHE-G and PAM50 scores
###########################################################

pattern <- "_full_cohort_predictions.csv"

compiled <- data.frame(matrix(nrow=0, ncol=0))
files <- list.files(dir, pattern = pattern)

for (file in files) {
    df <- read.csv(paste0(dir, "/", file))
    scores <- pam50_scores[match(df$sample, rownames(pam50_scores)),]
    arche <- toupper(sub(pattern, "", file))

    toCorr <- cbind(df, scores)
    tobind <- data.frame(
        Basal = cor.test(toCorr$y_pred, toCorr$Basal)$estimate,
        Her2 = cor.test(toCorr$y_pred, toCorr$Her2)$estimate,
        LumA = cor.test(toCorr$y_pred, toCorr$LumA)$estimate,
        LumB = cor.test(toCorr$y_pred, toCorr$LumB)$estimate,
        Normal = cor.test(toCorr$y_pred, toCorr$Normal)$estimate
    )
    rownames(tobind) <- arche
    compiled <- rbind(compiled, tobind)
}

###########################################################
# Plot compiled correlations
###########################################################

toPlot <- t(compiled)

# make colour palette
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
filename <- "data/results/figures/1-Signatures/ARCHE-G_subtype_heatmap.png"
png(filename, width = 3, height = 2, res = 600, units = "in")
ht
dev.off()

###########################################################
# Assign ARCHE-G
###########################################################

pattern <- "_full_cohort_predictions.csv"

compiled <- data.frame(matrix(nrow=0, ncol=0))
files <- list.files(dir, pattern = pattern)

for (file in files) {
    df <- read.csv(paste0(dir, "/", file))
    df$ARCHE <- toupper(sub(pattern, "", file))
    compiled <- rbind(compiled, df)
}

archeG <- data.frame(Sample = unique(compiled$sample))
archeG$ARCHEG <- archeG$Score <- NA

for (sample in archeG$Sample) {
    subset <- compiled[compiled$sample == sample,]
    top <- subset[subset$y_pred == max(subset$y_pred),]
    archeG$ARCHEG[archeG$Sample == sample] <- top$ARCHE
    archeG$Score[archeG$Sample == sample] <- top$y_pred
}
archeG$PAM50 <- pam50_scores$assigned[match(archeG$Sample, rownames(pam50_scores))]
saveRDS(archeG, file = "data/procdata/TCGA/ARCHEG_full_cohort.rds")

###########################################################
# Plot ARCHE-G PAM50 subtype counts
###########################################################

toPlot <- t(as.matrix(table(archeG$ARCHEG, archeG$PAM50)))

cols <- colorRampPalette(c("#EBEBEB", "#75BFB2", "#008A8C"))(9)
col_fun <- colorRamp2(seq(0,275,length.out = 9),cols)

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
filename <- "data/results/figures/1-Signatures/ARCHE-G_subtype_count_heatmap.png"
png(filename, width = 3, height = 2, res = 600, units = "in")
ht
dev.off()

###########################################################
# Plot survival plots
###########################################################

HR_res <- plot_survival(archeG, pheno, 5)