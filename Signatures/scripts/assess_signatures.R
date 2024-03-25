setwd("C:/Users/julia/Documents/BCaATAC")

library(umap)
library(ggplot2)
library(reshape2)
library(dplyr)

set.seed(101)

# set up palette for plotting
pal <- c("Signature1" = "#046C9A", "Signature2" = "#BBADB9", "Signature3" = "#7294D4", 
        "Signature4" = "#E8E1D9", "Signature5" = "#AFC5D8", "Signature6" = "#DF9C93")


# load in TCGA binary matrix
TCGA <- read.table("Signatures/data/BCa_binary.2.matrix", header = T)
nmf_scores <- read.table("Signatures/results/heatmap_rank6.png.order.matrix", header = T)

# extract consensus sequences
consensus <- TCGA[,1:3]
TCGA <- TCGA[,-c(1:3)]

# extract signatures
nmf_scores$Signature <- names(pal)
nmf_scores <- melt(nmf_scores)
nmf_scores <- as.data.frame(nmf_scores %>% group_by(variable) %>% filter(value == max(value)) %>% ungroup())

# transform matrix
TCGA <- t(TCGA)
sample_names <- rownames(TCGA)

# create UMAP projection
umap_df <- umap(TCGA)
umap_df <- as.data.frame(umap_df$layout)

# fromating umap df
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Signature <- nmf_scores[match(rownames(umap_df), nmf_scores$variable),]$Signature
umap_df <- rbind(umap_df, data.frame(UMAP1 = runif(20, min = -2, max = 2), UMAP2 = runif(20, min = -2, max = 2), Signature = ""))
#umap_df$sample <- sample_names
umap_df$type <- c(rep("Tumour", 75), rep("Cell Line", 10), rep("PDX", 10))
umap_df$type <- factor(umap_df$type, levels = c("Tumour", "Cell Line", "PDX"))

# plot umap
png("Signatures/results/signature_umap.png", width = 5, height = 4, res = 600, units = "in")
ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, shape = type)) + 
    geom_point(data = umap_df[umap_df$type == "Tumour",], aes(fill = Signature), size = 4, alpha = 0.8) +
    geom_point(data = umap_df[umap_df$type == "Cell Line",], fill = "#415A77", size = 2) +
    geom_point(data = umap_df[umap_df$type == "PDX",], fill = "#1B263B", size = 2) +
    scale_shape_manual(values = c(22, 21, 23)) +
    scale_fill_manual(values = c(pal)) +
    guides(shape = guide_legend(title = "Type", ncol = 1), fill = guide_legend(override.aes=list(shape = 22))) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
    labs(x = "UMAP1", y = "UMAP2", shape = "Type", fill = "Signature")
dev.off()
