setwd("C:/Users/julia/Documents/BCaATAC")

library(umap)
library(ggplot2)
library(reshape2)
library(dplyr)

set.seed(101)

# set up palette for plotting
pal <- c("Signature1" = "#046C9A", "Signature2" = "#BBADB9", "Signature3" = "#7294D4", 
        "Signature4" = "#E8E1D9", "Signature5" = "#AFC5D8", "Signature6" = "#DF9C93")


# load in TCGA cell line binary matrix
mat <- read.table("Signatures/data/BM_tcgabca_bcacelllines.tsv", header = T)

# load in TCGA binary matrix
nmf_scores <- read.table("Signatures/results/data/ATAC_heatmap_rank6.png.order.matrix", header = T)

# extract consensus sequences
consensus <- mat[,1:3]
mat <- mat[,-c(1:3)]

## TODO: Remove Lucie duplicates
# load in BCa cell lines
#samples <- read.csv("DrugResponse/data/cl_samples.csv")
#samples$file <- gsub("^X", "", gsub("\\.(?!$)", "-", samples$file, perl = TRUE))
#samples <- samples[match(umap_df$sample[76:138], samples$file),]
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

# extract signatures
nmf_scores$Signature <- paste0("Signature", 1:6)
nmf_scores <- melt(nmf_scores)
nmf_scores <- as.data.frame(nmf_scores %>% group_by(variable) %>% filter(value == max(value)) %>% ungroup())
tcga_signatures <- nmf_scores[match(colnames(mat)[1:75], nmf_scores$variable),]$Signature

# transform matrix
mat <- t(mat)
mat[is.na(mat)] <- 0
sample_names <- rownames(mat)

# create UMAP projection
umap_df <- as.data.frame(umap(mat)$layout)
save(umap_df, file = "Signatures/results/data/umap.RData")

# formating umap df
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Signature <- c(tcga_signatures, rep("", 63))
umap_df$type <- c(rep("Tumour", 75), rep("Cell Line", 63))
umap_df$type <- factor(umap_df$type, levels = c("Tumour", "Cell Line"))
umap_df$sample <- sample_names

save(umap_df, file = "Signatures/results/umap_df.RData")

# plot umap by signature
png("Signatures/results/figures/ATAC_umap_signature.png", width = 5, height = 4, res = 600, units = "in")
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




##### label by subtype

# set up palette for plotting
subtype_pal <- c("Basal" = "#AF4C5B","Her2" = "#EED4D3", "LumA" = "#B3B4D0", "LumB" = "#363E62", "Normal" = "#6365AF", "Not Available" = "#eFeBF7")


# load in TCGA annotation
meta <- read.csv("Signatures/data/TCGA_subtype_label.csv")
tumours <- gsub("\\.", "-", gsub("X", "", rownames(umap_df)[1:75]))
meta <- meta[match(tumours, meta$File.Name),]

# load in cell line annotation
cells_meta <- read.csv("MetaData/Lupien/BCa_samples.csv")
# from https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2911-z/tables/3, set 600MPE to LuminalA
cells_meta[which(cells_meta$subtype == "cell_line"),]$subtype <- "LumA"

# keep only cell lines being used
samples <- read.csv("DrugResponse/data/cl_samples.csv")
samples$match <- gsub("^X", "", gsub("\\.", "-", gsub("\\.$", "", samples$file)))
cells_meta <- cells_meta[cells_meta$sample %in% samples$match,]
cells <- rownames(umap_df)[76:138]
samples <- samples[match(cells, samples$file),]

# get subtypes
samples$subtype <- cells_meta[match(samples$match, cells_meta$sample),]$subtype
umap_df$subtype <- c(meta$bca_subtype, samples$subtype)


# plot umap by subtype
png("Signatures/results/figures/ATAC_umap_subtype.png", width = 5, height = 4, res = 600, units = "in")
ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, shape = type, fill = subtype)) + geom_point(size = 3) +
    scale_shape_manual(values = c(22, 21)) +
    scale_fill_manual(values = c(subtype_pal)) +
    guides(shape = guide_legend(title = "Type", ncol = 1), fill = guide_legend(override.aes=list(shape = 22))) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
    labs(x = "UMAP1", y = "UMAP2", shape = "Type", fill = "Subtype")
dev.off()





# Umap of just tumours
mat <- fread("Signatures/data/BCa_binary.2.matrix", header = T)

# extract consensus sequences and sample names
consensus <- mat[,1:3]
mat <- mat[,-c(1:3)]

# transform matrix
mat <- t(mat)
mat[is.na(mat)] <- 0

# create UMAP projection
umap_df <- as.data.frame(umap(mat)$layout)
save(umap_df, file = "Signatures/results/data/umap_tumours.RData")

load("Signatures/results/data/umap_tumours.RData")
colnames(umap_df) <- c("UMAP1", "UMAP2")

# set up palette for plotting
pal <- c("Signature1" = "#046C9A", "Signature2" = "#BBADB9", "Signature3" = "#7294D4", 
        "Signature4" = "#E8E1D9", "Signature5" = "#AFC5D8", "Signature6" = "#DF9C93")
subtype_pal <- c("Basal" = "#AF4C5B","Her2" = "#EED4D3", "LumA" = "#B3B4D0", "LumB" = "#363E62", "Normal" = "#6365AF", "Not Available" = "#eFeBF7")


# load in TCGA binary matrix
nmf_scores <- read.table("Signatures/results/data/ATAC_heatmap_rank6.png.order.matrix", header = T)

# extract signatures
nmf_scores$Signature <- paste0("Signature", 1:6)
nmf_scores <- melt(nmf_scores)
nmf_scores <- as.data.frame(nmf_scores %>% group_by(variable) %>% filter(value == max(value)) %>% ungroup())

# load in TCGA annotation
meta <- read.csv("Signatures/data/TCGA_subtype_label.csv")

# save signature and subtype annotation
umap_df$signature <- nmf_scores[match(rownames(umap_df), gsub("\\.", "-", gsub("X", "", nmf_scores$variable))),]$Signature
umap_df$subtype <- meta[match(rownames(umap_df), meta$File.Name),]$bca_subtype

# plot umap by subtype
png("Signatures/results/figures/ATAC_tumours_umap_subtype.png", width = 5, height = 4, res = 600, units = "in")
ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, fill = subtype)) + geom_point(size = 4, shape = 22) +
    scale_fill_manual(values = c(subtype_pal, "Not Available" = "white")) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
    labs(x = "UMAP1", y = "UMAP2", fill = "Subtype")
dev.off()

png("Signatures/results/figures/ATAC_tumours_umap_signature.png", width = 5.007, height = 4, res = 600, units = "in")
ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, fill = signature)) + geom_point(size = 4, shape = 22) +
    scale_fill_manual(values = c(pal)) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
    labs(x = "UMAP1", y = "UMAP2", fill = "Signature")
dev.off()
