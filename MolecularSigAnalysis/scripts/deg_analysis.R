setwd("Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(edgeR)
  library(openxlsx)
  library(ggplot2)
  library(ggrepel)
})


###########################################################
# Load in data
###########################################################

# load in counts matrix
counts <- fread("Signatures/data/TCGA_BRCA_gene_counts.matrix")

# load in metadata
meta <- read.csv("MolecularSigAnalysis/data/TCGA_sourcefiles.csv")

###########################################################
# Format counts matrix
###########################################################

samples <- counts$V1
counts$V1 <- NULL
counts <- t(counts)
colnames(counts) <- samples


###########################################################
# Designate groups
###########################################################

group <- meta$Signature[match(colnames(counts), meta$Sample.Name)]
group <- factor(group, levels = paste0("Signature", 1:6))


###########################################################
# Normalization and filtering
###########################################################

# create DGEList object
dge <- DGEList(counts=counts, group = group)

# filter genes with low counts
keep <- filterByExpr(dge, group = group)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# TMM normalization
dge <- calcNormFactors(dge)


###########################################################
# Differential expression
###########################################################

# function to perform DEG analysis per signature
DEG_sig <- function(signature) {
  
  # establish grouping of signatures
  if (signature > 1) {
    levels = paste0("Signature", c(signature:6, 1:(signature-1)) )
    sig_group = factor(group, levels = levels)
  } else {
    sig_group = group
  }
  design <- model.matrix(~sig_group)
  
  # DEG
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit, trend=TRUE)
  res <- topTable(fit, n=Inf, coef=ncol(design))
  
  # FDR 
  res$FDR <- p.adjust(res$P.Value, method = "BH", n = length(res$P.Value))
  
  return(res)
}

sig1 <- DEG_sig(1)
sig2 <- DEG_sig(2)
sig3 <- DEG_sig(3)
sig4 <- DEG_sig(4)
sig5 <- DEG_sig(5)
sig6 <- DEG_sig(6)



###########################################################
# Save top 15 genes
###########################################################

# function to extract top 15 genes
top15 <- function(sig_df) {
  top_genes <- sig_df[sig_df$FDR < 0.05,]
  top_genes <- top_genes[order(abs(top_genes$logFC), decreasing = T),]
  top_genes <- top_genes[1:15,]
  return(top_genes)
}

topSig1 <- top15(sig1)
topSig2 <- top15(sig2)
topSig3 <- top15(sig3)
topSig4 <- top15(sig4)
topSig5 <- top15(sig5)
topSig6 <- top15(sig6)

###########################################################
# Write dataframe results
###########################################################

path = "MolecularSigAnalysis/results/data/deg_"

write.xlsx(sig1, file = paste0(path, "signature1", ".xlsx"))
write.xlsx(sig2, file = paste0(path, "signature2", ".xlsx"))
write.xlsx(sig3, file = paste0(path, "signature3", ".xlsx"))
write.xlsx(sig4, file = paste0(path, "signature4", ".xlsx"))
write.xlsx(sig5, file = paste0(path, "signature5", ".xlsx"))
write.xlsx(sig6, file = paste0(path, "signature6", ".xlsx"))

write.csv(topSig1, file = paste0(path, "top_signature1", ".csv"))
write.csv(topSig2, file = paste0(path, "top_signature2", ".csv"))
write.csv(topSig3, file = paste0(path, "top_signature3", ".csv"))
write.csv(topSig4, file = paste0(path, "top_signature4", ".csv"))
write.csv(topSig5, file = paste0(path, "top_signature5", ".csv"))
write.csv(topSig6, file = paste0(path, "top_signature6", ".csv"))

###########################################################
# Volcano plots
###########################################################

# set palette for plotting
pal = c("Upregulated" = "#93E1D8", "Downregulated" = "#AA4465")

# function to plot DEG and label top 15 genes
volcano <- function(sig_df, label, top_genes) {
  
  top_genes$dir <- ifelse(top_genes$logFC > 0, "Upregulated", "Downregulated")
  top_genes$dir <- factor(top_genes$dir, levels = c("Upregulated", "Downregulated"))
  
  ggplot() + 
    geom_point(data = sig_df, aes(x = logFC, y = FDR), color = "gray") +
    geom_point(data = top_genes, aes(x = logFC, y = FDR, color = dir)) +
    geom_text_repel(data = top_genes, max.overlaps = 100, size = 2.5,
                    aes(x = logFC, y = FDR, label = rownames(top_genes))) +
    scale_color_manual("", values = pal) +
    ylim(1:-0.15) + xlim(-1.55,1.55) + 
    theme_classic() + theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(label)
}

png("MolecularSigAnalysis/results/figures/volcano_sig1.png", width = 8, height = 6, res = 600, units = "in")
volcano(sig1, "Signature1", topSig1)
dev.off()

png("MolecularSigAnalysis/results/figures/volcano_sig2.png", width = 8, height = 6, res = 600, units = "in")
volcano(sig2, "Signature2", topSig2)
dev.off()

png("MolecularSigAnalysis/results/figures/volcano_sig3.png", width = 8, height = 6, res = 600, units = "in")
volcano(sig3, "Signature3", topSig3)
dev.off()

png("MolecularSigAnalysis/results/figures/volcano_sig4.png", width = 8, height = 6, res = 600, units = "in")
volcano(sig4, "Signature4", topSig4)
dev.off()

png("MolecularSigAnalysis/results/figures/volcano_sig5.png", width = 8, height = 6, res = 600, units = "in")
volcano(sig5, "Signature5", topSig5)
dev.off()

png("MolecularSigAnalysis/results/figures/volcano_sig6.png", width = 8, height = 6, res = 600, units = "in")
volcano(sig6, "Signature6", topSig6)
dev.off()


###########################################################
# Plot top 15 genes
###########################################################

# function to plot top 15 genes
plot_top15 <- function(top15, label) {
  
  top15 <- top15[order(top15$logFC, decreasing = F),]
  top15$rank <- 1:15
  top15$Gene <- rownames(top15)
  top15$Gene <- factor(top15$Gene, levels = top15$Gene)
  
  p <- ggplot(top15) + geom_tile(aes(x = 1, y = Gene, fill = logFC), color = "gray") +
    scale_fill_gradient2(high = "#93E1D8", low = "#AA4465", mid = "white", midpoint = 0) +
    geom_text(aes(x = 1, y = Gene, label = Gene), color = "black", size = 3) +
    theme_void() + 
    theme(axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
    labs(x = "", y = "Gene") + theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(label)
  return(p)
}

png("MolecularSigAnalysis/results/figures/topSig2.png", width = 2, height = 4, res = 600, units = "in")
plot_top15(topSig2, "Signature2")
dev.off()

png("MolecularSigAnalysis/results/figures/topSig3.png", width = 2, height = 4, res = 600, units = "in")
plot_top15(topSig3, "Signature3")
dev.off()

png("MolecularSigAnalysis/results/figures/topSig4.png", width = 2, height = 4, res = 600, units = "in")
plot_top15(topSig4, "Signature4")
dev.off()

png("MolecularSigAnalysis/results/figures/topSig5.png", width = 2, height = 4, res = 600, units = "in")
plot_top15(topSig5, "Signature5")
dev.off()

png("MolecularSigAnalysis/results/figures/topSig6.png", width = 2, height = 4, res = 600, units = "in")
plot_top15(topSig6, "Signature6")
dev.off()
