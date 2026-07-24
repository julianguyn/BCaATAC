# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(matrixStats)
    library(patchwork)
    library(ggnewscale)
})

source("utils/get_data.R")
source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# read in meta data file
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")

# load in tumour gene counts matrix
t_counts <- get_tcga_rna()

###########################################################
# Load in gene signatures
###########################################################

# read in lineage genes
genes <- read.csv("data/rawdata/lineages/43018_2024_773_MOESM2_ESM.xlsx - S5 Supplementary Table 5.csv", skip = 1)
genes$gene <- sub("\\..*", "", genes$gene)
#Basal_myoepithelial      Luminal_mature  Luminal_progenitor 
#               1213                1225                 927

# subset to top genes
genes <- genes[abs(genes$avg_log2FC) > 1,]
#Basal_myoepithelial      Luminal_mature  Luminal_progenitor 
#                138                 152                  94

genes$direction <- ifelse(genes$avg_log2FC > 0, "Upregulated", "Downregulated")
#table(genes$direction, genes$cluster)
#                Basal_myoepithelial Luminal_mature Luminal_progenitor
#  Downregulated                  64             87                 45
#  Upregulated                    74             65                 49


# get gene signatures
bm_up <- genes$gene[genes$cluster == "Basal_myoepithelial" & genes$direction == "Upregulated"]
bm_dn <- genes$gene[genes$cluster == "Basal_myoepithelial" & genes$direction == "Downregulated"]
lm_up <- genes$gene[genes$cluster == "Luminal_mature" & genes$direction == "Upregulated"]
lm_dn <- genes$gene[genes$cluster == "Luminal_mature" & genes$direction == "Downregulated"]
lp_up <- genes$gene[genes$cluster == "Luminal_progenitor" & genes$direction == "Upregulated"]
lp_dn <- genes$gene[genes$cluster == "Luminal_progenitor" & genes$direction == "Downregulated"]

###########################################################
# Helper function to score gene signatures
###########################################################

score_signature <- function(rna, up_genes, down_genes, label) {
  
  # missing genes
  up_present   <- intersect(up_genes, rownames(rna))
  down_present <- intersect(down_genes, rownames(rna))
  
  missing_up   <- setdiff(up_genes, rownames(rna))
  missing_down <- setdiff(down_genes, rownames(rna))
  
  if (length(missing_up) > 0) {
    warning(sprintf("%d up-genes not found in expression matrix: %s",
                     length(missing_up), paste(missing_up, collapse = ", ")))
  }
  if (length(missing_down) > 0) {
    warning(sprintf("%d down-genes not found in expression matrix: %s",
                     length(missing_down), paste(missing_down, collapse = ", ")))
  }
  if (length(up_present) == 0 | length(down_present) == 0) {
    stop("No overlapping up- or down-genes found in expression matrix.")
  }
  
  # z-score genes across samples (row-wise)
  gene_means <- rowMeans(rna, na.rm = TRUE)
  gene_sds <- rowSds(rna, na.rm = TRUE)
  gene_sds[gene_sds == 0] <- NA
  z_mat <- (rna - gene_means) / gene_sds
  
  # mean z-score per sample
  up_scores <- colMeans(z_mat[up_present,,drop = FALSE], na.rm = TRUE)
  down_scores <- colMeans(z_mat[down_present,,drop = FALSE], na.rm = TRUE)
  
  res <- data.frame(
    sample = colnames(rna),
    up_score = up_scores,
    down_score = down_scores,
    sig_score = up_scores - down_scores,
    label = label,
    row.names  = NULL
  )
  return(res)
}

###########################################################
# Score signatures
###########################################################

score <- rbind(
    score_signature(t_counts, bm_up, bm_dn, "Basal_myoepithelial"),
    score_signature(t_counts, lm_up, lm_dn, "Luminal_mature"),
    score_signature(t_counts, lp_up, lp_dn, "Luminal_progenitor")
)
score$ARCHE <- meta$ARCHE[match(score$sample, meta$Sample.Name)]
score$Subtype <- meta$Subtype[match(score$sample, meta$Sample.Name)]

###########################################################
# Plot scores
###########################################################

plot_lineage_score <- function(score, label) {

    toPlot <- score[score$label == label,]
    p <- ggplot(toPlot, aes(x = ARCHE, y = sig_score)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
        geom_boxplot(aes(fill = ARCHE), outlier.shape = NA) + 
        scale_fill_manual(values = ARCHE_pal) +
        new_scale_fill() +
        geom_jitter(aes(fill = Subtype), 
                    shape = 21, colour = "black", size = 2, width = 0.3
                    ) +
        scale_fill_manual("Subtype", values = subtype_pal) +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.title.x = element_blank()
        ) +
        labs(y = paste(label, "\nscore"))
    return(p)
}

p1 <- plot_lineage_score(score, "Basal_myoepithelial")
p2 <- plot_lineage_score(score, "Luminal_mature")
p3 <- plot_lineage_score(score, "Luminal_progenitor")

p <- p1 / p2 / p3
ggsave("data/results/figures/1-Signatures/ARCHE_lineages.png", p, width = 5, height = 7)
