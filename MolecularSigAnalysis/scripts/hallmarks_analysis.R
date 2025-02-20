setwd("Documents/BCaATAC")


# load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
  library(RColorBrewer)
})


###########################################################
# Load in data
###########################################################

# load in counts matrix
counts <- fread("Signatures/data/TCGA_BRCA_gene_counts.matrix") |>
  as.data.frame()

# load hallmarks gene set
msigh <- msigdbr(species = "human", category = "H")

# load in matrix file from NMF
mat <- read.table("ATAC_heatmap_rank6.png.order.matrix", header = T) |>
  as.data.frame()

# read in meta data file
meta <- read.csv("TCGA_sourcefiles.csv")
meta$ATAC.Seq.File.Name <- gsub("-", "\\.", meta$ATAC.Seq.File.Name)


###########################################################
# Format counts matrix and hallmarks gene sets
###########################################################

# format counts matrix
rownames(counts) <- counts$V1
counts$V1 <- NULL

# format hallmarks gene sets
hallmarks <- split(x = msigh$gene_symbol, f = msigh$gs_name)


###########################################################
# GSEA on hallmark gene sets
###########################################################

# initiate dataframe to store results
results <- data.frame(hallmarks = unique(msigh$gs_name))

# extract hallmark enrichment scores for each sample
for (i in seq_along(rownames(counts))) {
  print(i)
  ranks = counts[i,,drop = F] |> as.numeric()
  names(ranks) <- colnames(counts)
  fgseaRes <- fgsea(pathways = hallmarks, stats = ranks)
  results <- cbind(results, fgseaRes$ES)
}
colnames(results) <- c("hallmarks", rownames(counts))


###########################################################
# Save results
###########################################################

write.csv(results, file = "MolecularSigAnalysis/results/data/hallmarks_ES.csv", quote = F, row.names = F)


###########################################################
# Format matrices for correlation
###########################################################

# format ATAC signature matrix
rownames(mat) <- paste0("Signature", 1:6)

# rename samples
colnames(mat) <- meta$Sample.Name[match(gsub("X", "", colnames(mat)), meta$ATAC.Seq.File.Name)]

# format hallmarks matrix
rownames(results) <- results$hallmarks
results$hallmarks <- NULL


###########################################################
# Keep common samples
###########################################################

common <- intersect(colnames(mat), colnames(results))
results <- results[,match(common, colnames(results))]
mat <- mat[,match(common, colnames(mat))]


###########################################################
# Correlate each pair of signatures
###########################################################

# initiate dataframe to hold results
corr <- data.frame(matrix(nrow=0, ncol=3))
colnames(corr) <- c("ATAC.Sig", "Hallmark", "Corr")

for (i in seq_along(rownames(mat))) {
  
  atac.sig <- rownames(mat)[i]
  
  for (j in seq_along(rownames(results))) {
    
    hallmark <- rownames(results)[j]
    s <- cor(as.numeric(mat[i,]), as.numeric(results[j,]), method = "spearman")
    corr <- rbind(corr, data.frame(ATAC.Sig = atac.sig, Hallmark = hallmark, Corr = s))
  }
}


###########################################################
# Plot heatmap of correlations
###########################################################

# format for plotting
corr$ATAC.Sig <- factor(corr$ATAC.Sig, levels = c(paste0("Signature", 1:6)))

# plot heatmap
png("MolecularSigAnalysis/results/figures/hallmarks.png", width = 6, height = 8, res = 600, units = "in")
ggplot(corr) + geom_tile(aes(x = Hallmark, y = ATAC.Sig, fill = Corr), color = "gray") +
  scale_fill_gradient2("Spearman\nCorrelation",
                       low = "#AF4C5B", high = "#046C9A", mid = "white",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_void() + 
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y = element_text(size = 8, hjust = 0.95),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
  labs(x = "Hallmark Gene Sets", y = "\nATAC Signatures") + coord_flip()
dev.off()