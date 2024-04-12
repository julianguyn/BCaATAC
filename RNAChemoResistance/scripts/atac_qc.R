# attempt batch correction on atac
setwd("C:/Users/julia/Documents/BCaATAC")

suppressMessages(library(sva))
library(ggplot2)

# load in loci counts from Jocelyn
df <- read.csv("RNAChemoResistance/data/rawCounts.csv")
colnames(df) <- gsub("\\.", "-", colnames(df))

# format counts df
rownames(df) <- df$X
df <- df[,-c(1)]

# read in metadata
meta <- read.csv("RNAChemoResistance/data/config.csv")
meta <- meta[match(colnames(df), meta$SampleID),]

# remove loci with low counts (SKIP BC NO LOW COUNTS)
#exp_mean <- apply(df,1,mean)
#mat <- as.matrix(df[-which(rownames(df) %in% names(exp_mean[which(unname(exp_mean) == 0)])),])

# remove outlier
outlier <- "MDAMB436_sens_3"
df <- df[,-which(colnames(df) == outlier)]
meta <- meta[-which(meta$SampleID == outlier),]

mat <- df

# adjust for batch
count_adj <- ComBat(dat=mat, batch = meta$Factor, par.prior=TRUE, prior.plots=FALSE)
#count_adj_seq <- ComBat_seq(mat, meta$Factor, group=NULL)
write.table(count_adj, file = "RNAChemoResistance/results/data/atac_countsadj.tsv", quote = F, col.names = TRUE, sep = "\t")


# view with PCA
pca_org <- as.data.frame(prcomp(t(df))$x)
pca_adj <- as.data.frame(prcomp(t(count_adj))$x)

# add in labels for plotting PCA scores
pca_org$batch <- meta$Factor
pca_adj$batch <- meta$Factor

pca_org$label <- meta$Condition
pca_adj$label <- meta$Condition

# plot PCA to show batches by cell
png("RNAChemoResistance/results/figures/atac_pca_org.png", width = 8, height = 6, res = 600, units = "in")
ggplot(data = pca_org, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = label, shape = batch), size = 3) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(title = "PCA Original")
dev.off()

png("RNAChemoResistance/results/figures/atac_pca_adj.png", width = 8, height = 6, res = 600, units = "in")
ggplot(data = pca_adj, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = label, shape = batch), size = 3) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 15)) +
  labs(title = "PCA Adjusted")
dev.off()
