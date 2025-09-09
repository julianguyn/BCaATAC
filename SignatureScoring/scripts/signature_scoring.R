setwd("/home/bioinf/bhklab/julia/projects/ATACseq")

library(data.table)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

# set the number of signatures
num_signatures <- 6

# set colours for plotting
subtype_pal <- c("Basal" = "#AF4C5B","Her2" = "#EED4D3", "LumA" = "#B3B4D0", "LumB" = "#363E62", "Normal" = "#6365AF", "Not Available" = "#eFeBF7")


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


## unknown files
cells <- as.data.frame(fread("Signature_Scoring/data/bcaCelllines.matrix", header = T))
pdxs <- as.data.frame(fread("Signature_Scoring/data/pdx.matrix", header = T))

#################################
######## Tumour Heatmaps ########
#################################

#########################
### ATAC-Seq Profiles ###

# set up palette for plotting
pal <- c("Signature1" = "#046C9A", "Signature2" = "#BBADB9", "Signature3" = "#7294D4", 
        "Signature4" = "#E8E1D9", "Signature5" = "#AFC5D8", "Signature6" = "#DF9C93",
        "Signature7" = "#4E4C67", "Signature8" = "#B4869F", "Signature9" = "#A6B1E1",
        "Signature10" = "#985F6F")
subtype_pal <- c("Basal" = "#394032","Her2" = "#A6A57A", "LumA" = "#8C5E58", 
                 "LumB" = "#5A352A", "Normal" = "#8F8073", "Not Available" = "#eFeBF7")

# read in meta data
meta <- read.csv("Signature_Scoring/data/TCGA_subtype_label.csv")

# load in matrix file from NMF
scores <- read.table("Signature_Scoring/data/ATAC_heatmap_rank6.png.order.matrix", header = T)
scores$Signature <- paste0("Signature", 1:6)
scores <- melt(scores)
scores$variable <- gsub("^X", "", scores$variable)

# load in binary matrix
mat <- as.data.frame(fread("Signature_Scoring/data/BCa_binary.2.matrix", header = T))

rownames(mat) <- paste(mat$V1, mat$V2, mat$V3, sep = ":")
mat <- mat[,-c(1:3)]

# compute correlations
corr_df <- cor(mat, use = "complete.obs")


# set colours for plotting
score_pal = colorRamp2(seq(min(corr_df), max(corr_df), length = 3), c("#B1EBEE", "#348AA7", "#513B56"))

# get signature assignment
signatures <- c()
for (sample in gsub("-", "\\.", colnames(corr_df))) {
    tmp <- scores[scores$variable == sample,]
    signatures <- c(signatures, as.character(tmp[which.max(tmp$value),]$Signature))
}

# set subtype and clustering order annotation
col_fun <- colorRamp2(c(0, 30, 50), c("#CEF0F2", "#42B6A3", "#446D73"))
ha = HeatmapAnnotation(Subtype = meta[match(colnames(corr_df), meta$File.Name),]$bca_subtype, Signature = signatures,
                       col = list(Subtype = subtype_pal, Signature = pal))

# plot heatmap
ht <- Heatmap(corr_df, name = "Correlation", col = score_pal,
        column_title = "Tumour Sample", column_title_side = "bottom", show_column_names = FALSE, show_row_names = FALSE,
        top_annotation = ha)
ht <- draw(ht)

png("Signature_Scoring/results/figures/ATAC_tumours_heatmap.png", width = 9, height = 8, res = 600, units = "in")
ht
dev.off()



####################################
######## Cell Line Heatmaps ########
####################################

# load in metadata
meta <- fread("Signature_Scoring/data/bcacells_annotation.tsv")
samples <- read.csv("Signature_Scoring/data/cl_samples.csv")

# remove duplicates
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

# save subtype annotation
samples$subtype <- meta[match(samples$sample, meta$sample),]$subtype

#########################
### ATAC-Seq Profiles ###

# load in counts matrix
mat <- as.data.frame(fread("Signature_Scoring/data/bcaCelllines.matrix", header = T))

rownames(mat) <- paste(mat$chr, mat$start, mat$end, sep = ":")
mat <- mat[,-c(1:3)]
colnames(mat) <- gsub("-", "\\.", colnames(mat))

# keep only needed samples
mat <- mat[,colnames(mat) %in% gsub("\\.$", "", samples$file)]
colnames(mat) <- samples[match(colnames(mat), gsub("\\.$", "", samples$file)),]$sample

# compute deviations
mat2 <- as.data.frame(t(compute_dev(as.data.frame(t(mat)))))

# compute correlations
corr_df <- cor(mat, use = "complete.obs")


# set colours for plotting
score_pal = colorRamp2(seq(-0.5, 1, length = 3), c("#B1EBEE", "#348AA7", "#513B56"))

# get clustering order
ht <- Heatmap(corr_df, name = "Correlation", col = score_pal)
ht <- draw(ht)
atac_order <- row_order(ht)
order_df <- data.frame(sample = colnames(corr_df)[atac_order], order = 1:49)
save(order_df, file = "SignatureScoring/results/cell_line_order.RData")

# set subtype and clustering order annotation
col_fun <- colorRamp2(c(0, 30, 50), c("#CEF0F2", "#42B6A3", "#446D73"))
ha = HeatmapAnnotation(Subtype = samples[match(colnames(corr_df), samples$sample),]$subtype, Order = order_df[match(colnames(corr_df), order_df$sample),]$order,
                       col = list(Subtype = subtype_pal, Order = col_fun))

# plot heatmap
ht <- Heatmap(corr_df, name = "Correlation", col = score_pal,
        column_title = "Cell Lines", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
        top_annotation = ha)
ht <- draw(ht)

png("Signature_Scoring/results/figures/ATAC_cells_heatmap.png", width = 9, height = 8, res = 600, units = "in")
ht
dev.off()




########################
### RNA-Seq Profiles ###
########################

# load in gene matrix
mat <- as.data.frame(fread("Signature_Scoring/data/bcacells_gene_counts.matrix", header = T))

rownames(mat) <- mat[,c(1)]
mat <- mat[,-c(1)]
mat <- as.data.frame(t(mat))

# match to sample sheet
samples$subtype <- meta[match(samples$sample, meta$sample),]$subtype
samples$order <- order_df[match(samples$sample, order_df$sample),]$order
samples$match <- gsub("-", "", samples$sample)

# compute correlations
corr_df <- cor(mat)

# set colours for plotting
score_pal = colorRamp2(seq(0, 1, length = 3), c("#B1EBEE", "#348AA7", "#513B56"))

# set subtype and clustering order annotation
col_fun <- colorRamp2(c(0, 30, 50), c("#CEF0F2", "#42B6A3", "#446D73"))
ha = HeatmapAnnotation(Subtype = samples[match(colnames(corr_df), samples$match),]$subtype, Order = samples[match(colnames(corr_df), samples$match),]$order,
                       col = list(Subtype = subtype_pal, Order = col_fun))

# plot heatmap
ht <- Heatmap(corr_df, name = "Correlation", col = score_pal,
        column_title = "Cell Lines", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 9),
        top_annotation = ha)
ht <- draw(ht)

png("Signature_Scoring/results/figures/RNA_cells_heatmap.png", width = 9, height = 8, res = 600, units = "in")
ht
dev.off()



# Quick umap of RNA-Seq
library(umap)
library(ggplot2)

umap_df <- as.data.frame(umap(t(mat))$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

png("Signature_Scoring/results/figures/RNA_cells_umap.png", width = 5, height = 4, res = 600, units = "in")
ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2)) + geom_point(size = 2.5) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
    labs(x = "UMAP1", y = "UMAP2")
dev.off()

ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, shape = type, fill = subtype)) + geom_point(size = 2.5) +
    scale_shape_manual(values = c(22, 21)) +
    scale_fill_manual(values = c(subtype_pal)) +
    guides(shape = guide_legend(title = "Type", ncol = 1), fill = guide_legend(override.aes=list(shape = 22))) +
    theme_classic() + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.key.size = unit(0.7, 'cm')) +
    labs(x = "UMAP1", y = "UMAP2", shape = "Type", fill = "Subtype")



###########################
### ATAC-Seq Signatures ###
###########################

# load in ChromVar output
scores <- read.table("Signature_Scoring/data/bca_sign.Zscore.txt")

# keep only samples of interest
scores <- scores[,which(colnames(scores) %in% samples$file)]
colnames(scores) <- samples[match(colnames(scores), samples$file),]$sample
rownames(scores) <- paste0("Signature", 1:6)

# normalize scores
scores <- 1 - compute_dev(scores)
scores <- (scores * -1) + 1
saveRDS(scores, file = "Signature_Scoring/results/data/cl_signaturescores.RDS")

# remove two cell lines (CIHR)
scores <- scores[,-which(colnames(scores) %in% c("CAL-51", "Hs 578T"))]
rownames(scores) <- paste0("CAS", 1:6)

# set colours for plotting
score_pal = colorRamp2(seq(min(scores), max(scores), length = 3), c("#C3BFCC", "#F8F1F8", "#077293"))

# get subtype
subtype_pal <- c("Basal" = "#AF4C5B","Her2" = "#EED4D3", "LumA" = "#B3B4D0", "LumB" = "#363E62", "Normal" = "#6365AF", "Not Available" = "#eFeBF7")

ha = HeatmapAnnotation(Subtype = meta[match(colnames(scores), meta$sample),]$subtype, 
                       col = list(Subtype = subtype_pal))

png("SignatureScoring/results/figures/ATAC_cells_signatures_heatmap.png", width = 8, height = 4, res = 600, units = "in")
Heatmap(scores, cluster_rows = FALSE, name = "CAS\nExpression\nScore", col = score_pal,
        column_title = "Cell Lines", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 10),
        top_annotation = ha)
dev.off()


####################################
# Subtype per CAS for CIHR 
####################################

# Subtype per CAS for CIHR
meta <- fread("DataExploration/data/bcacells_annotation.tsv")

scores$CAS <- rownames(scores)
test <- melt(scores) |> suppressWarnings()
test$Subtype <- meta[match(test$variable, meta$sample),]$subtype

#ggplot(test, aes(x = Subtype, y = value)) + geom_boxplot() + facet_wrap(~ CAS)

png("SignatureScoring/results/figures/subtype_per_CAS.png", width = 5, height = 5.5, res = 600, units = "in")
ggplot(test, aes(x = value, y = Subtype, fill = Subtype)) + geom_boxplot() + 
        facet_grid(CAS ~., switch = "y") +
        scale_fill_manual(values = subtype_pal) +
        theme_classic() +
        geom_vline(xintercept = 0.5, linetype = "dashed") +
        theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
              legend.key.size = unit(0.7, 'cm'),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), 
              panel.spacing = unit(0, "pt")) +
        labs(y = "", x = "CAS Expression Score")
dev.off()


##########################
### RNA-Seq Signatures ###
##########################

# load in counts 
count_df <- as.data.frame(fread("Signature_Scoring/data/bcacells_gene_counts.matrix"))

# format count matrix
rownames(count_df) <- count_df$Gene
count_df <- count_df[,-c(1)]


# load in gene set (signatures)
gene_set <- as.data.frame(fread("Signatures/results/RNA_NMF_output/rank6Mixture.csv"))
gene_anno <- fread("Signatures/results/RNA_NMF_output/gene_annotation.txt", header = F)
colnames(gene_set) <- gene_anno$V1
rownames(gene_set) <- paste0("Signature", 1:num_signatures)

# keep only intersected genes
intersected <- intersect(colnames(gene_set), colnames(count_df)) #16783
gene_set <- gene_set[,colnames(gene_set) %in% intersected]
count_df <- count_df[,colnames(count_df) %in% intersected]
count_df <- count_df[,match(colnames(gene_set), colnames(count_df))]

# normalize data
gene_set <- compute_dev(gene_set)
count_df <- compute_dev(count_df)

# compute 1 - Pearson's correlation
corr_df <- 1 - cor(t(gene_set), t(count_df), use = "complete.obs")

# set colours for plotting
score_pal = colorRamp2(seq(min(corr_df), max(corr_df), length = 3), c("#148BAF", "white", "#A06B9A"))

# get subtype
meta <- meta[which(meta$match %in% colnames(corr_df)),]
ha = HeatmapAnnotation(Subtype = meta[match(colnames(corr_df), meta$match),]$subtype, 
                       col = subtype_pal)


png("Signature_Scoring/results/figures/RNA_cells_signatures_heatmap.png", width = 8, height = 4, res = 600, units = "in")
Heatmap(corr_df, cluster_rows = FALSE, name = "Signature\nScore", col = score_pal,
        column_title = "Cell Lines", column_title_side = "bottom", column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 10),
        top_annotation = ha)
dev.off()


