# load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
})

source("utils/plots/deg_analysis.R")
source("utils/get_data.R")
source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# load in counts matrix
counts <- get_tcga_rna()

# save samples w rna data (for heatmap)
have_rna <- gsub("\\.", "-", colnames(counts))
saveRDS(have_rna, file = "data/procdata/TCGA/samples_w_rna.RDS")

# load in metadata file
meta <- read.csv("data/rawdata/tcga/TCGA_sourcefiles.csv")
meta <- meta[match(colnames(counts), meta$Sample.Name),]

###########################################################
# Un-log2 counts
###########################################################

# counts matrix was log2 normalized, undo that
counts <- round(2^counts - 1)

###########################################################
# Function for DEG
###########################################################

run_DEG <- function(counts, meta, arche = NA, subtype = NA) {

  # create colData
  colData <- data.frame(
    SampleID = meta$Sample.Name,
    ARCHE = ifelse(meta$ARCHE == arche, arche, "Other"),
    Subtype = ifelse(meta$Subtype == subtype, subtype, "Other")
  )
  
  if (!is.na(arche)) {  # if running ARCHEx vs all
    colData$ARCHE <- factor(colData$ARCHE, levels = c("Other", arche))
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design= ~ ARCHE) |>
                                suppressMessages()
    label <- paste0("ARCHE_", arche, "_vs_Other")
  } else {            # if running Subtypex vs all
    colData$Subtype <- factor(colData$Subtype, levels = c("Other", subtype))
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design= ~ Subtype) |>
                                suppressMessages()
    label <- paste0("Subtype_", subtype, "_vs_Other")
  }

  # DEG
  dds <- DESeq(dds) |> suppressMessages()
  res <- results(dds, name=label)
  write.csv(res, file = paste0("data/results/data/2-MolecularSigAnalysis/DEG/", label, ".csv"), quote = F, row.names = T)
  return(res)
}

###########################################################
# Run DEG of each ARCHE vs all
###########################################################

a1 <- run_DEG(counts, meta, arche = "ARCHE1")
a2 <- run_DEG(counts, meta, arche = "ARCHE2")
a3 <- run_DEG(counts, meta, arche = "ARCHE3")
a4 <- run_DEG(counts, meta, arche = "ARCHE4")
a5 <- run_DEG(counts, meta, arche = "ARCHE5")
a6 <- run_DEG(counts, meta, arche = "ARCHE6")

# ---------------------------------------------------------
# Load DEGs
# ---------------------------------------------------------

dir <- "data/results/data/2-MolecularSigAnalysis/DEG/"

a1 <- read.csv(paste0(dir, "ARCHE_ARCHE1_vs_Other.csv"))
a2 <- read.csv(paste0(dir, "ARCHE_ARCHE2_vs_Other.csv"))
a3 <- read.csv(paste0(dir, "ARCHE_ARCHE3_vs_Other.csv"))
a4 <- read.csv(paste0(dir, "ARCHE_ARCHE4_vs_Other.csv"))
a5 <- read.csv(paste0(dir, "ARCHE_ARCHE5_vs_Other.csv"))
a6 <- read.csv(paste0(dir, "ARCHE_ARCHE6_vs_Other.csv"))


###########################################################
# Subset and format DEGs
###########################################################

subset_sig <- function(deg) {
  sig <- deg[which(abs(deg$log2FoldChange) > 2 & deg$padj < 0.05),]
  return(sig)
}

a1 <- subset_sig(a1)
a2 <- subset_sig(a2)
a3 <- subset_sig(a3)
a4 <- subset_sig(a4)
a5 <- subset_sig(a5)
a6 <- subset_sig(a6)

sig_genes <- unique(c(a1$X, a2$X, a3$X, a4$X, a5$X, a6$X))

###########################################################
# Format data for heatmap of upregulated genes
###########################################################

# assign each upregulated gene to ARCHE where it is most upregulated
degs <- list(ARCHE1 = a1, ARCHE2 = a2, ARCHE3 = a3, ARCHE4 = a4, ARCHE5 = a5, ARCHE6 = a6)
up <- do.call(rbind, lapply(names(degs), function(a) {
  d <- degs[[a]][degs[[a]]$log2FoldChange > 0,]
  data.frame(gene = d$X, ARCHE = a, lfc = d$log2FoldChange)
}))
up <- up[order(-up$lfc),]
up <- up[!duplicated(up$gene),]

toPlot <- counts[up$gene,]
gene_group <- factor(up$ARCHE, levels = names(ARCHE_pal))

# order samples by ARCHE
meta <- meta[order(meta$ARCHE), ]
toPlot <- toPlot[, match(meta$Sample.Name, colnames(toPlot))]

# log2 and z-score each gene across samples
toPlot <- t(scale(t(log2(toPlot + 1))))

# order genes by ARCHE group
gene_order <- unlist(lapply(levels(gene_group), function(a) {
  idx <- which(gene_group == a)
  idx[hclust(dist(toPlot[idx,]))$order]
}))
toPlot <- toPlot[gene_order,]
gene_group <- gene_group[gene_order]

###########################################################
# Make heatmap
###########################################################

# make colour palette
cols <- colorRampPalette(c("#1F1636", "#30244E", "#50466B", "#9CC6C0", "#C2DFDB"))(9)
col_fun <- colorRamp2(
    seq(-2,2,
        length.out = 9),
    cols
)

ha <- HeatmapAnnotation(
  'ARCHE' = meta$ARCHE,
  col = list('ARCHE' = ARCHE_pal),
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(3, "mm")
)

# genes to label on heatmap
shortlist <- list(
  ARCHE1 = c("NPY1R", "RERG", "CITED1"),
  ARCHE2 = c("FOXC1", "ELF5", "VGLL1", "SOX10", "KRT5", "KRT17", "CRYAB", "ACTR3B"),
  ARCHE3 = c("ERBB2", "C17orf37", "GSDMB", "PPP1R1B", "S100A8", "DMBT1"),
  ARCHE4 = c("IGF1R", "TRH"),
  ARCHE5 = c("SPIB", "CD274", "LAG3", "IDO1", "GZMB", "IFNG", "CXCL10", "VCAM1"),
  ARCHE6 = c("FOXJ1", "CDC20B", "IGFALS", "CXCL14")
)
label_idx <- which(rownames(toPlot) %in% unlist(shortlist))

row_ha <- rowAnnotation(
  'Genes' = anno_mark(
    at = label_idx,
    labels = rownames(toPlot)[label_idx],
    side = "left",
    labels_gp = gpar(fontsize = 9),
    link_width = unit(3, "mm"),
    link_gp = gpar(lwd = 0.75)
  ),
  'Gene ARCHE' = gene_group,
  col = list('Gene ARCHE' = ARCHE_pal),
  simple_anno_size = unit(3, "mm"),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht <- Heatmap(
  toPlot,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  left_annotation = row_ha
)

png("data/results/figures/2-MolecularSigAnalysis/deg/ARCHE_sig_genes_heatmap.png", width = 5, height = 7, res = 600, units = "in")
draw(ht, merge_legend = TRUE)

# diagonal boxes
row_ends <- cumsum(table(gene_group)) / nrow(toPlot)
col_ends <- cumsum(table(factor(meta$ARCHE, levels = levels(gene_group)))) / ncol(toPlot)
decorate_heatmap_body("Z-score", {
  for (i in seq_along(row_ends)) {
    x0 <- c(0, col_ends)[i]; y0 <- c(0, row_ends)[i]
    grid.rect(
      x = x0, y = 1 - y0, width = col_ends[i] - x0, height = row_ends[i] - y0,
      just = c("left", "top"), gp = gpar(fill = NA, col = "white", lwd = 1.5)
    )
  }
})

dev.off()


###########################################################
# Volcano plots
###########################################################

plot_volcano(a1, "ARCHE1")
plot_volcano(a2, "ARCHE2")
plot_volcano(a3, "ARCHE3")
plot_volcano(a4, "ARCHE4")
plot_volcano(a5, "ARCHE5")
plot_volcano(a6, "ARCHE6")

plot_volcano(sb, "Basal")
plot_volcano(sh, "Her2")
plot_volcano(sla, "LumA")
plot_volcano(slb, "LumB")
#plot_volcano(sn, "Normal") # only 4 sig genes, too small sample size

###########################################################
# Volcano of ARCHE5 TFs
###########################################################

tfs <- c("SPIB", "IRF3", "IRF8", "SPI1")
a2$Label <- ifelse(rownames(a2) %in% tfs, "Immune TF", "Other")
a2$Gene <- rownames(a2)

p <- ggplot(a2, aes(x = log2FoldChange, y = -log(padj))) +
  geom_point(color = "gray") +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "#989898") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#989898") +
  geom_point(data = a2[a2$Label == "Immune TF", ], 
             aes(x = log2FoldChange, y = -log(padj)), color = "red") +
  geom_text_repel(data = a2[a2$Label == "Immune TF", ],
                   aes(x = log2FoldChange, y = -log(padj), label = Gene),
                   color = "black", size = 3, max.overlaps = Inf) +
  theme_bw()

ggsave("data/results/figures/2-MolecularSigAnalysis/deg/ARCHE2_TF.png", p, width = 4.5, height = 4)

###########################################################
# Plot MYC volcano plots
###########################################################

plot_volcano_MYC(a2, "ARCHE2")
plot_volcano_MYC(a5, "ARCHE5")
plot_volcano_MYC(sb, "Basal")

###########################################################
# MYC expression analysis across all tumours
###########################################################

myc <- as.data.frame(counts[rownames(counts) == "MYC",]) 
myc$Signature <- meta$Signature[match(rownames(myc), meta$Sample.Name)]
myc$Subtype <- meta$Subtype[match(rownames(myc), meta$Sample.Name)]
myc$Sample.Name <- rownames(myc)
colnames(myc)[1] <- "MYC"

# plot MYC expression
plot_MYCexp(myc, "Signature", "all")
plot_MYCexp(myc, "Subtype", "all")

###########################################################
# ARCHE 2 vs ARCHE 5 MYC expression
###########################################################

# myc expression 
myc <- myc[which(myc$Subtype == "Basal" & myc$Signature != "ARCHE3"),]
plot_MYCexp(myc, "Signature", "basal")

# DEG
myc_counts <- counts[,match(rownames(myc), colnames(counts))]
a25 <- run_DEG(myc_counts, myc, arche = "ARCHE2")
plot_volcano(a25, "Basal ARCHE2")
plot_volcano_MYC(a25, "Basal ARCHE2")

###########################################################
# Label DMR genes in DEG volcano plots
###########################################################

# load in genes from DMRs (ARCHE2 and 3)
load("data/results/data/2-MolecularSigAnalysis/DMR_genes.RData")

plot_volcano_select_genes(a2, arche2_genes, "ARCHE2", subset = TRUE)
plot_volcano_select_genes(a2, overlap, "MYC_ARCHE2")
plot_volcano_select_genes(a3, arche3_genes, "ARCHE3", subset = TRUE)