# compile tumour fractions for REFLECT from ichorCNA

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
})

source("utils/palettes.R")

set.seed(101)

###########################################################
# Load in data
###########################################################

# load in metadata
meta <- read_excel("data/rawdata/cfDNA/REFLECT-unfiltered/REFLECT-6B_metadata.xlsx", sheet = 1) |> as.data.frame()
meta$id_6b <- gsub("_", "-", gsub("_LB.*", "", meta$library_id))
meta$Subtype_final[meta$Subtype_final == "ER+/HER2-"] <- "ER+"
meta$Subtype_final[meta$Subtype_final == "TBNC"] <- "TNBC"
meta$Subtype_final[meta$Subtype_final == "HER2+"] <- "HER2"

###########################################################
# Get tumour fraction estimates + other stats
###########################################################

dir <- "data/rawdata/cfdna-ichorCNA/TF/"
results <- data.frame(matrix(nrow=7, ncol=0))

# loop through params.txt
for (file in list.files(dir)) {
    params <- readLines(paste0(dir, file))
    split_params <- strsplit(params[6:11], "\t")
    vals <- c(params[5], sapply(split_params, `[`, 2))
    results[[vals[1]]] <- vals
}
results <- t(results) |> as.data.frame()
colnames(results) <- c("Sample", "Gender", "Tumour_Fraction", "Ploidy", "Subclone_Fraction", "Fraction_Genome_Subclonal", "Fraction_CNA_Subclonal")

# get subtype
results$subtype <- meta$Subtype_final[match(results$Sample, meta$id_6b)]

# format results
results$Tumour_Fraction <- as.numeric(results$Tumour_Fraction)
results$Ploidy <- as.numeric(results$Ploidy)
results$Subclone_Fraction <- as.numeric(results$Subclone_Fraction)
results$Fraction_Genome_Subclonal <- as.numeric(results$Fraction_Genome_Subclonal)
results$Fraction_CNA_Subclonal <- as.numeric(results$Fraction_CNA_Subclonal)

# format results for plotting
results <- results[order(results$Tumour_Fraction, decreasing = TRUE),]
results$Sample <- factor(results$Sample, levels = results$Sample)
results$REF052 <- ifelse(results$Sample == "REFLECT-0052-03", "REF052", "Other")

# save results
write.csv(results, file = "data/results/data/5-cfDNA/REFLECT/tumour_fractions.csv", quote = FALSE, row.names = FALSE)

###########################################################
# Plot counts
###########################################################

# get numbers
n5 <- nrow(results[results$Tumour_Fraction > 0.05,])
n10 <- nrow(results[results$Tumour_Fraction > 0.1,])

p <- ggplot(results, aes(x = Sample, y = Tumour_Fraction, fill = REF052)) +
    geom_bar(stat = "identity", color = "black") +
    geom_hline(yintercept = c(0.05, 0.1), linetype = "dashed") +
    annotate("text", x = Inf, y = 0.05, label = paste("TF > 5%:", n5), 
             hjust = 1.1, vjust = -0.5, size = 3) +
    annotate("text", x = Inf, y = 0.1, label = paste("TF > 10%:", n10), 
             hjust = 1.1, vjust = -0.5, size = 3) +
    annotate("text", x = Inf, y = 0.6, label = paste("Number of Samples", nrow(results)), 
             hjust = 1.1, vjust = -0.5, size = 3) +
    scale_fill_manual(values = REF052_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm')
    )
filename <- paste0("data/results/figures/Misc/ichorCNA/REFLECT_TF_scores.png")
png(filename, width=8, height=5, units='in', res = 600, pointsize=80)
p
dev.off()

###########################################################
# Plot heatmap components
###########################################################

results$Sample <- factor(results$Sample, levels = rev(results$Sample))

# counts plot
p <- ggplot(results, aes(x = Tumour_Fraction, y = Sample, fill = REF052)) +
    geom_bar(stat = "identity", color = "black") +
    geom_vline(xintercept = c(0.05, 0.1), linetype = "dashed") +
    scale_fill_manual(values = REF052_pal) +
    theme_minimal() +
    theme(
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, 'cm')
    ) 
filename <- paste0("data/results/figures/Misc/ichorCNA/REFLECT_heatmap2.png")
png(filename, width=5, height=8, units='in', res = 600, pointsize=80)
p
dev.off()

# heatmap
toPlot <- results[,c(grep("Fraction", colnames(results)))]

# set colours for plotting
score_pal = colorRamp2(seq(0, 1, length = 3), c(binary_pal2[2], "#F8F1F8", binary_pal2[1]))

ha <- rowAnnotation(
  Subtype = results$subtype[match(rownames(toPlot), rownames(results))],
  REF052 = results$REF052[match(rownames(toPlot), rownames(results))],
  col = list(Subtype = cfDNA_subtype_pal, REF052 = REF052_pal))

filename <- paste0("data/results/figures/Misc/ichorCNA/REFLECT_heatmap1.png")
png(filename, width = 5, height = 10, res = 600, units = "in")
Heatmap(toPlot, cluster_rows = FALSE, name = "TF Binding\nScore", col = score_pal,
  column_title = "Samples", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 10), right_annotation = ha)
dev.off()
