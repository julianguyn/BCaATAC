# Script to just explore the Griffin output

library(data.table)
library(reshape2)
library(matrixStats)

# read in files
CAMA1_C6M1 <- fread("CAMA1_C6M1.GC_corrected.coverage.tsv") |> as.data.frame()
CAMA1_C6M4 <- fread("CAMA1_C6M4.GC_corrected.coverage.tsv") |> as.data.frame()
CAMA1_PAR_Rep1 <- fread("CAMA1_PAR_Rep1.GC_corrected.coverage.tsv") |> as.data.frame()
CAMA1_PAR_Rep2 <- fread("CAMA1_PAR_Rep2.GC_corrected.coverage.tsv") |> as.data.frame()
T47D_PAR_Rep1 <- fread("T47D_PAR_Rep1.GC_corrected.coverage.tsv") |> as.data.frame()
T47D_PAR_Rep2 <- fread("T47D_PAR_Rep2.GC_corrected.coverage.tsv") |> as.data.frame()
ZR751_PAR_Rep1 <- fread("ZR751_PAR_Rep1.GC_corrected.coverage.tsv") |> as.data.frame()
ZR751_PAR_Rep2 <- fread("ZR751_PAR_Rep2.GC_corrected.coverage.tsv") |> as.data.frame()

# read in data frame
# dataframe from 'Calculating coverage signature scores at top 1000 sig sites' from Nucleosome_Profile_Analysis.ipynb
cov <- rbind(c(0.688349, 0.696604, 0.731168, 0.702734, 0.728619, 0.747513, 0.623232, 0.575476), 
             c(1.003672, 0.977843, 0.997243, 1.014506, 1.007018, 1.006096, 0.905857, 1.074052),
             c(0.862610, 0.879886, 1.004805, 0.951594, 0.804883, 0.807255, 0.873123, 0.907878),
             c(0.940813, 0.928023, 0.973092, 0.949821, 1.002566, 0.971887, 0.870868, 0.929523),
             c(0.973293, 0.957655, 1.022560, 0.971874, 0.986281, 0.981508, 1.016437, 0.934534),
             c(0.859333, 0.827164, 0.863944, 0.886805, 0.841030, 0.836796, 0.777311, 0.896960)) |>
    as.data.frame()
colnames(cov) <- c("CAMA1_Rep1", "CAMA1_Rep2", "T47D_Rep1", "T47D_Rep2",
                   "ZR751_Rep1", "ZR751_Rep2", "CAMA1_C6M1", "CAMA1_C6M4")
cov$Signature <- paste0("Signature", 1:6)

# format dataframe for plotting
toPlot <- melt(cov)
toPlot$Signature <- factor(toPlot$Signature, levels = paste0("Signature", 6:1))

# plot heatmap
png("cfDNA/results/figures/griffin_coverage.png", width=125, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = variable, y = Signature, fill = value)) + 
    geom_tile(color = "white") + theme_void() +
    geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
    scale_fill_gradientn(colors = c("#077293", "#77A1AE", "#F8F1F8")) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.5, hjust = 1),
        axis.title.y = element_text(angle = 90, size = 12),
        strip.text.x = element_text(),
        legend.key.size = unit(0.7, 'cm')
    )  +
    labs(x = "Sample", y = "Signature\n", fill = "Mean\nCoverage")
dev.off()

# plot without mice samples
toPlot <- melt(cov[,-which(colnames(cov) %in% c("CAMA1_C6M1", "CAMA1_C6M4"))])
toPlot$Signature <- factor(toPlot$Signature, levels = paste0("Signature", 6:1))

# plot heatmap
png("cfDNA/results/figures/griffin_coverage_cellsonly.png", width=100, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = variable, y = Signature, fill = value)) + 
    geom_tile(color = "white") + theme_void() +
    geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
    scale_fill_gradientn(colors = c("#077293", "#77A1AE", "#F8F1F8")) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.5, hjust = 1),
        axis.title.y = element_text(angle = 90, size = 12),
        strip.text.x = element_text(),
        legend.key.size = unit(0.7, 'cm')
    )  +
    labs(x = "Sample", y = "Signature\n", fill = "Mean\nCoverage")
dev.off()



# replot ATAC signature scores
scores <- readRDS("SignatureScoring/results/cl_signaturescores.RDS")

# keep only cell lines of interest
cl_scores <- data.frame(
    'CAMA-1.r1' = scores$CAMA-1,
    'CAMA-1.r2' = scores$CAMA-1,
    'T-47D.r1' = scores$'T-47D',
    'T-47D.r2' = scores$'T-47D',
    'ZR-75-1.r1' = scores$ZR-75-1,
    'ZR-75-1.r2' = scores$ZR-75-1
)

# functon to compute deviations from mean to normalize data
compute_dev <- function(df) {
    means <- colMeans(df)
    stdev <- colSds(as.matrix(df))

    df[is.na(df)] <- 0

    norm_df <- as.data.frame(sapply(1:ncol(df), function(i) { (df[[i]] - means[i]) / stdev[i] }))
    colnames(norm_df) <- colnames(df)
    rownames(norm_df) <- rownames(df)

    return(norm_df)
}

mat <- compute_dev(cl_scores)
mat$Signature <- paste0("Signature", 1:6)

# format dataframe for plotting
toPlot <- melt(mat)
toPlot$Signature <- factor(toPlot$Signature, levels = paste0("Signature", 6:1))

# plot heatmap
png("cfDNA/results/figures/atac_sig_score.png", width=100, height=150, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = variable, y = Signature, fill = value)) + 
    geom_tile(color = "white") + theme_void() +
    geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
    scale_fill_gradientn(colors = c("#C3BFCC", "#F8F1F8", "#077293"), limits = c(-2, 2)) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.5, hjust = 1),
        axis.title.y = element_text(angle = 90, size = 12),
        strip.text.x = element_text(),
        legend.key.size = unit(0.7, 'cm')
    )  +
    labs(x = "Sample", y = "Signature\n", fill = "Signature\nScore")
dev.off()


### ===== Correlation between CAMA1 cell line and mice ===== ### 

# subset to CAMA1 cell line and mice
cama <- cov[,grep("CAMA", colnames(cov))]

# create correlation matrix
cor_mat <- cor(cama)

# plot correlation
toPlot <- melt(cor_mat)

png("cfDNA/results/figures/corr_cama.png", width=90, height=80, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile(color = "white") + theme_void() +
    geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
    scale_fill_gradientn(colors = c("#C3BFCC", "#F8F1F8", "#077293"), limits = c(-1, 1)) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.5, hjust = 1),
        axis.title.y = element_text(angle = 90, size = 12),
        strip.text.x = element_text(),
        legend.key.size = unit(0.7, 'cm')
    )  +
    labs(x = "", y = "", fill = "Pearson")
dev.off()

### ===== Correlation between cfDNA and ATAC ===== ### 

# remove Signature column
cov <- cov[,-which(colnames(cov) == "Signature")]
mat <- mat[,-which(colnames(mat) == "Signature")]

# remove CAMA-1 mice samples
cov <- cov[,-grep("C6M", colnames(cov))]

# inverse cfDNA signature scores
cov <- -cov

# spearman correlations per sample
sample_res <- data.frame(matrix(nrow = ncol(cov), ncol = 2))
colnames(sample_res) <- c("Sample", "Spearman")

for (i in 1:ncol(cov)) {
    sample = colnames(cov)[i]

    x <- cor.test(cov[,i], mat[,i], method = "spearman") |> suppressWarnings()
    sample_res$Sample[i] <- sample
    sample_res$Spearman[i] <- x$estimate
}

# plot correlation
toPlot <- melt(sample_res)
toPlot$Sample <- factor(toPlot$Sample, levels = rev(toPlot$Sample))

png("cfDNA/results/figures/corr_sample.png", width=60, height=80, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = 1, y = Sample, fill = value)) + 
    geom_tile(color = "white") + theme_void() +
    geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
    scale_fill_gradientn(colors = c("#BC4749", "#F8F1F8", "#077293"), limits = c(-1, 1)) +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.5, hjust = 1),
        axis.title.y = element_text(angle = 90, size = 12),
        strip.text.x = element_text(),
        legend.key.size = unit(0.7, 'cm')
    )  +
    labs(x = "Spearman", y = "", fill = "Spearman")
dev.off()

# spearman correlations per signature
signature_res <- data.frame(matrix(nrow = nrow(cov), ncol = 2))
colnames(signature_res) <- c("Signature", "Spearman")

for (i in 1:nrow(cov)) {

    x <- cor.test(as.numeric(cov[i,]), as.numeric(mat[i,]), method = "spearman") |> suppressWarnings()
    signature_res$Signature[i] <- paste0("Signature", i)
    signature_res$Spearman[i] <- x$estimate
}

# plot correlation
toPlot <- melt(signature_res)
toPlot$Signature <- factor(toPlot$Signature, levels = rev(toPlot$Signature))

png("cfDNA/results/figures/corr_signature.png", width=55, height=80, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = 1, y = Signature, fill = value)) + 
    geom_tile(color = "white") + theme_void() +
    geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
    scale_fill_gradientn(colors = c("#BC4749", "#F8F1F8", "#077293"), limits = c(-1, 1)) +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(vjust = 0.5, hjust = 1),
        axis.title.y = element_text(angle = 90, size = 12),
        strip.text.x = element_text(),
        legend.key.size = unit(0.7, 'cm')
    )  +
    labs(x = "Spearman", y = "", fill = "Spearman")
dev.off()