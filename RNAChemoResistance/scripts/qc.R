setwd("C:/Users/julia/Documents/BCaATAC")

suppressMessages(library(sva))
library(ggplot2)

# load in counts and PCA from Jocelyn
load("RNAChemoResistance/data/counts.RData")
load("RNAChemoResistance/data/PCA.RData")

# read in metadata
meta <- read.csv("RNAChemoResistance/data/metadata.csv")
pca_scores$batch <- meta[match(rownames(pca_scores), paste0(meta$Sampleid, meta$Rep)),]$Batch

# format counts matrix and remove genes with low counts
rownames(count_df) <- count_df$Geneid
count_df <- count_df[,-c(1)]
exp_mean <- apply(count_df,1,mean)
mat <- as.matrix(count_df[-which(rownames(count_df) %in% names(exp_mean[which(unname(exp_mean) == 0)])),])

# adjust for batch
count_adj <- ComBat_seq(mat, batch=pca_scores$batch, group=NULL)
write.table(count_adj, file = "RNAChemoResistance/results/data/adjcounts.tsv", quote = F, col.names = TRUE, sep = "\t")

# view with PCA
pca_adj <- as.data.frame(prcomp(t(count_adj))$x)

# add in labels for plotting PCA scores
pca_adj$cell <- ifelse(rownames(pca_adj) %in% rownames(pca_adj)[grep("HS578T", rownames(pca_adj))], "HS578T", "MDAMB436")
pca_adj$label <- gsub("\\d", "", gsub("MDAMB436", "", gsub("HS578T", "", rownames(pca_adj))))
pca_adj$label <- ifelse(pca_adj$label %in% c("A", "B", "C", ""), paste0("sens",pca_adj$label), pca_adj$label)
pca_adj$drug <- gsub("[ABC]", "", pca_adj$label)
pca_adj$batch <- meta[match(rownames(pca_adj), paste0(meta$Sampleid, meta$Rep)),]$Batch

# plot PCA to show batches by cell
png("RNAChemoResistance/results/figures/pca_adj_cell.png", width = 8, height = 6, res = 600, units = "in")
ggplot(data = pca_adj, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = cell, shape = batch), size = 3) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(title = "PCA Adjusted")
dev.off()

# plot PCA to show batches by label
png("RNAChemoResistance/results/figures/pca_adj_lab.png", width = 8, height = 6, res = 600, units = "in")
ggplot(data = pca_adj, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = label, shape = batch), size = 3) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(title = "PCA Adjusted")
dev.off()



##################################################################################
# split by cell type
MDA_counts <- count_df[,grep("MDA", colnames(count_df))]
HS_counts <- count_df[,grep("HS578T", colnames(count_df))]
gene_names <- rownames(count_df)

# perform PCA
MDA_pca <- as.data.frame(prcomp(t(MDA_counts))$x)
HS_pca <- as.data.frame(prcomp(t(HS_counts))$x)

# add labels for plotting PCA scores
MDA_pca$label <- gsub("\\d", "", gsub("MDAMB436", "", rownames(MDA_pca)))
MDA_pca$label <- ifelse(MDA_pca$label %in% c("A", "B", "C"), paste0("sens",MDA_pca$label), MDA_pca$label)
MDA_pca$drug <- gsub("[ABC]", "", MDA_pca$label)
MDA_pca$batch <- meta[match(rownames(MDA_pca), paste0(meta$Sampleid, meta$Rep)),]$Batch

png("RNAChemoResistance/results/figures/pca_MDA.png", width = 8, height = 6, res = 600, units = "in")
ggplot(data = MDA_pca, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = label, shape = batch), size = 3) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(title = "MDAMB436 PCA Original")
dev.off()

###############################
#### Adjusting for Batches ####
###############################

# format counts matrix and remove genes with low counts
rownames(MDA_counts) <- rownames(count_df)
exp_mean <- apply(MDA_counts,1,mean)
mat <- as.matrix(MDA_counts[-which(rownames(MDA_counts) %in% names(exp_mean[which(unname(exp_mean) == 0)])),])

# adjust for batch
MDA_counts_adj <- ComBat_seq(mat, batch=MDA_pca$batch, group=NULL)
write.table(MDA_counts_adj, file = "RNAChemoResistance/results/data/adjcounts.tsv", quote = F, col.names = TRUE, sep = "\t")

# view with PCA
MDA_pca_adj <- as.data.frame(prcomp(t(MDA_counts_adj))$x)

# add in labels for plotting PCA scores
MDA_pca_adj$label <- gsub("\\d", "", gsub("MDAMB436", "", rownames(MDA_pca_adj)))
MDA_pca_adj$label <- ifelse(MDA_pca_adj$label %in% c("A", "B", "C"), paste0("sens",MDA_pca_adj$label), MDA_pca_adj$label)
MDA_pca_adj$drug <- gsub("[ABC]", "", MDA_pca_adj$label)
MDA_pca_adj$batch <- meta[match(rownames(MDA_pca_adj), paste0(meta$Sampleid, meta$Rep)),]$Batch

# plot PCA to show batches by cell
png("RNAChemoResistance/results/figures/pca_MDA_adj.png", width = 8, height = 6, res = 600, units = "in")
ggplot(data = MDA_pca_adj, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = label, shape = batch), size = 3) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(title = "MDAMB436 PCA Adjusted")
dev.off()


###########################################
#### Visualizating Counts Distribution ####
###########################################

library(ggplot2)
library(tidyr)

# plot original 
genes_all <- gather(as.data.frame(mat), key = "Sample", value = "Expression")
genes_all$Expression <- log2(genes_all$Expression + 1)

png("RNAChemoResistance/results/figures/MDA_density.png", width = 15, height = 6, res = 600, units = "in")
ggplot(genes_all, aes(x = Expression, color = Sample)) +
  geom_density() + theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(title = "Density Plot of Gene Expression by Sample", x = "Expression", y = "Density")
dev.off()

# plot adjusted for cell 
genes_all <- gather(as.data.frame(MDA_counts_adj), key = "Sample", value = "Expression")
genes_all$Expression <- log2(genes_all$Expression + 1)

png("RNAChemoResistance/results/figures/MDA_density_adj.png", width = 15, height = 6, res = 600, units = "in")
ggplot(genes_all, aes(x = Expression, color = Sample)) +
  geom_density() + theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(title = "Density Plot of Gene Expression by Sample", x = "Expression", y = "Density")
dev.off()