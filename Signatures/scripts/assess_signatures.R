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
nmf_scores <- read.table("Signatures/results/ATAC_heatmap_rank6.png.order.matrix", header = T)

# load in cell line binary matrix
cells <- read.table("Signatures/data/bcacell_lines.tsv", header = T)

# extract consensus sequences
consensus <- TCGA[,1:3]
cell_peaks <- cells[,1:3]

# find overlap in peak calls
consensus$merge <- paste0(consensus$V1, "-", consensus$V2, "-", consensus$V3)
cell_peaks$merge <- paste0(cell_peaks$V1, "-", cell_peaks$V2, "-", cell_peaks$V3)

consensus$cells <- 0

# subset by chromosome for quicker search
subset <- cell_peaks[cell_peaks$V1 == "chr1",]

for (i in 1:nrow(consensus)) {
    
    # get new subset if chr changes
    if (subset$V1[1] != consensus$V1[i]) {
        print(paste("---------starting", consensus$V1[i]))
        subset <- cell_peaks[cell_peaks$V1 == consensus$V1[i], ]
    }
    if (i %% 5000 == 0) {print(paste(i, "out of", 627305))} 

    if (consensus$merge[i] %in% subset$merge) {

        # find direct overlap
        consensus$cells[i] <- subset$merge[i] 

    } else {

        # find the closest peak
        subset$diff <- abs(subset$V2 - consensus$V2[i]) + abs(subset$V3 - consensus$V3[i])
        closest <- subset$merge[which.min(subset$diff)]
        consensus$cells[i] <- closest
    }
}
print("saving")
save(consensus, file = "Signatures/results/TCGA_cell_match.RData")

## TODO: Remove Lucie duplicates
# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")
samples <- samples[match(umap_df$sample[76:138], samples$file),]
# TEMP
# dup <- samples$sample[duplicated(samples$sample)]
# tmp <- samples[samples$sample %in% dup,]
# tmp <- tmp[tmp$seq == "Lucie",]
# samples$label <- ifelse(samples$file %in% tmp$file, "Lucie-toRemove", samples$seq) # from dups, keep only the ones by Nergiz
# umap_df$cell_label <- c(rep("Tumour", 75), samples$label)
# umap_df$cell_label <- factor(umap_df$cell_label, levels = c("Tumour", "Nergiz", "Lucie", "Lucie-toRemove"))
#png("Signatures/results/signature_umap_Lucie.png", width = 5, height = 4, res = 600, units = "in")
#ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, shape = cell_label)) + 
#     geom_point(data = umap_df[umap_df$type == "Tumour",], aes(fill = Signature), size = 4, alpha = 0.8) +
#     geom_point(data = umap_df[umap_df$type == "Cell Line",], fill = "#415A77", size = 2) +
#     #geom_point(data = umap_df[umap_df$type == "PDX",], fill = "#1B263B", size = 2) +
#     #scale_shape_manual(values = c(22, 21, 23)) +
#     scale_shape_manual(values = c(22, 23, 8, 1)) +
#     scale_fill_manual(values = c(pal)) +
#     guides(shape = guide_legend(title = "Type", ncol = 1), fill = guide_legend(override.aes=list(shape = 22))) +
#     theme_classic() + 
#     theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
#     labs(x = "UMAP1", y = "UMAP2", shape = "Type", fill = "Signature")
#dev.off()

# reorder cell peak calls to TCGA consensus sequence
cells_matched <- cells[match(consensus$cells, cell_peaks$merge),]
rownames(cells_matched) <- NULL

# combine TCGA tumour and cell binary files
TCGA <- TCGA[,-c(1:3)]
cells_matched <- cells_matched[,-c(1:3)]

combined <- cbind(TCGA, cells_matched)

# extract signatures
nmf_scores$Signature <- paste0("Signature", 1:6)
nmf_scores <- melt(nmf_scores)
nmf_scores <- as.data.frame(nmf_scores %>% group_by(variable) %>% filter(value == max(value)) %>% ungroup())
tcga_signatures <- nmf_scores[match(colnames(TCGA), nmf_scores$variable),]$Signature

# transform matrix
combined <- t(combined)
combined[is.na(combined)] <- 0
sample_names <- rownames(combined)

# create UMAP projection
umap_df <- umap(combined)
umap_df <- as.data.frame(umap_df$layout)

# formating umap df
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Signature <- c(tcga_signatures, rep("", ncol(cells_matched)))
umap_df$type <- c(rep("Tumour", ncol(TCGA)), rep("Cell Line", ncol(cells_matched)))
umap_df$type <- factor(umap_df$type, levels = c("Tumour", "Cell Line"))
umap_df$sample <- sample_names

save(umap_df, file = "Signatures/results/umap_df.RData")

# plot umap
png("Signatures/results/signature_umap.png", width = 5, height = 4, res = 600, units = "in")
ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, shape = type)) + 
    geom_point(data = umap_df[umap_df$type == "Tumour",], aes(fill = Signature), size = 4, alpha = 0.8) +
    geom_point(data = umap_df[umap_df$type == "Cell Line",], fill = "#415A77", size = 2) +
    #geom_point(data = umap_df[umap_df$type == "PDX",], fill = "#1B263B", size = 2) +
    #scale_shape_manual(values = c(22, 21, 23)) +
    scale_shape_manual(values = c(22, 21)) +
    scale_fill_manual(values = c(pal)) +
    guides(shape = guide_legend(title = "Type", ncol = 1), fill = guide_legend(override.aes=list(shape = 22))) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
    labs(x = "UMAP1", y = "UMAP2", shape = "Type", fill = "Signature")
dev.off()

