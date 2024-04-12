# Make UMAP of RNA-Seq gene counts between tumour and cell lines coloured by subtype

setwd("C:/Users/julia/Documents/BCaATAC")

library(umap)
library(ggplot2)
library(reshape2)
library(dplyr)
library(PharmacoGx)

set.seed(101)

# set up palette for plotting
subtype_pal <- c("Basal" = "#394032", "Her2" = "#A6A57A","LumA" = "#8C5E58","LumB" = "#5A352A","Normal" = "#8F8073", "#eFeBF7")

# load in TCGA gene counts matrix
TCGA <- read.table("Signatures/data/TCGA_BRCA_gene_counts.matrix", header = T)

# load in TCGA annotation
meta <- read.csv("Signatures/data/TCGA_subtype_label.csv")
meta <- unique(meta[,c(3:4)])
meta <- meta[match(rownames(TCGA), meta$sample_name),]

# load in cell line RNA seq counts
uhnbreast2 <- readRDS("DrugResponse/data/PharmacoSet.RDS")
cells <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@assays@data$expr
colnames(cells) <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@colData$sampleid
genes_meta <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@rowRanges
rownames(cells) <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@rowRanges$gene_name

# load in cell line annotation
cells_meta <- read.csv("MetaData/Lupien/BCa_samples.csv")
# from https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2911-z/tables/3, set 600MPE to LuminalA
cells_meta[which(cells_meta$subtype == "cell_line"),]$subtype <- "LumA"

# keep only cell lines being used
samples <- read.csv("DrugResponse/data/cl_samples.csv")
samples$file <- gsub("\\.$", "", samples$file)
samples$file <- gsub("\\.", "-", samples$file)
samples$match <- gsub("-", "", samples$sample)

# remove dups
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

cells_meta <- cells_meta[which(cells_meta$sample %in% samples$file),]
cells <- cells[,colnames(cells) %in% samples$match]

# format dataframe
cells <- as.data.frame(t(cells))

# save labels
tcga_labs <- rownames(TCGA)
cell_labs <- rownames(cells)

# keep common genes
keep <- intersect(colnames(TCGA), colnames(cells))
TCGA <- TCGA[,colnames(TCGA) %in% keep]
cells <- cells[,colnames(cells) %in% keep]
cells <- cells[,match(colnames(TCGA), colnames(cells))]

# combine dataframes
combined <- rbind(TCGA, cells)

# create UMAP projection
umap_df <- as.data.frame(umap(combined)$layout)

# formating umap df
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$type <- c(rep("Tumour", nrow(TCGA)), rep("Cell Line", nrow(cells)))
umap_df$type <- factor(umap_df$type, levels = c("Tumour", "Cell Line"))
umap_df$sample <- c(tcga_labs, cell_labs)

# get molecular subtype
cells_meta$match <- samples[match(cells_meta$sample, samples$file),]$match
umap_df$subtype <- c(meta[match(tcga_labs, meta$sample_nam),]$bca_subtype, cells_meta[match(cell_labs, cells_meta$match),]$subtype)

# plot umap
png("DataExploration/results/figures/rna_umap.png", width = 5, height = 4, res = 600, units = "in")
ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, shape = type, fill = subtype)) + geom_point(size = 2.5) +
    scale_shape_manual(values = c(22, 21)) +
    scale_fill_manual(values = c(subtype_pal)) +
    guides(shape = guide_legend(title = "Type", ncol = 1), fill = guide_legend(override.aes=list(shape = 22))) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
    labs(x = "UMAP1", y = "UMAP2", shape = "Type", fill = "Subtype")
dev.off()

