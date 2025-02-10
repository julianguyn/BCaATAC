setwd("Downloads/maf")

# load libraries
suppressPackageStartupMessages({
  library(maftools)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  library(ggh4x)
  #library(BSgenome.Hsapiens.UCSC.hg19)
})

###########################################################
# Load in data
###########################################################

# load in matrix file from NMF
mat <- read.table("../ATAC_heatmap_rank6.png.order.matrix", header = T)

# read in meta data file
meta <- read.csv("MolecularSigAnalysis/data/TCGA_sourcefiles.csv")
meta$ATAC.Seq.File.Name <- gsub("-", "\\.", meta$ATAC.Seq.File.Name)
meta$SNV.File.Name <- gsub("\\.gz", "", meta$SNV.File.Name)

###########################################################
# Get signature assignment
###########################################################

# get signature
mat$Signature <- paste0("Signature", 1:6)

# format for plotting
mat <- reshape2::melt(mat)

# save signature assignment
mat$signature_assign <- ""
for (sample in mat$variable) {
  tmp <- mat[mat$variable == sample,]
  mat[mat$variable == sample,]$signature_assign <- as.character(tmp[which.max(tmp$value),]$Signature)
}

# keep only variables of interest
mat <- mat[,colnames(mat) %in% c("variable", "signature_assign")] |> unique()

# save signature to metadata
mat$variable <- gsub("X", "", mat$variable)
meta$Signature <- mat$signature_assign[match(meta$ATAC.Seq.File.Name, mat$variable)]

# keep order of samples
mat$rank = 1:nrow(mat)
meta$rank <- mat$rank[match(meta$ATAC.Seq.File.Name, mat$variable)]


###########################################################
# MAF summaries all tumours
###########################################################

# load in all maf files
mafs <- merge_mafs(meta$SNV.File.Name |> na.omit())

png("MolecularSigAnalysis/results/figures/maf.png", width = 8, height = 6, res = 600, units = "in")
plotmafSummary(maf = mafs, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

###########################################################
# MAF summaries for each signature
###########################################################

# function to generate maf summary
signature_mafsummary <- function(signature) {
  files <- meta$SNV.File.Name[meta$Signature == signature] |> na.omit()
  mafs <- merge_mafs(files)
  p <- plotmafSummary(maf = mafs, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  return(p)
}

png("MolecularSigAnalysis/results/figures/sig1_maf.png", width = 8, height = 6, res = 600, units = "in")
signature_mafsummary("Signature1")
dev.off()

png("MolecularSigAnalysis/results/figures/sig2_maf.png", width = 8, height = 6, res = 600, units = "in")
signature_mafsummary("Signature2")
dev.off()

png("MolecularSigAnalysis/results/figures/sig3_maf.png", width = 8, height = 6, res = 600, units = "in")
signature_mafsummary("Signature3")
dev.off()

png("MolecularSigAnalysis/results/figures/sig4_maf.png", width = 8, height = 6, res = 600, units = "in")
signature_mafsummary("Signature4")
dev.off()

png("MolecularSigAnalysis/results/figures/sig5_maf.png", width = 8, height = 6, res = 600, units = "in")
signature_mafsummary("Signature5")
dev.off()

png("MolecularSigAnalysis/results/figures/sig6_maf.png", width = 8, height = 6, res = 600, units = "in")
signature_mafsummary("Signature6")
dev.off()

###########################################################
# Create mutations counts matrix
###########################################################

# create counts matrix
mat <- mutCountMatrix(
  mafs,
  includeSyn = FALSE,
  countOnly = NULL,
  removeNonMutated = TRUE
)


###########################################################
# Map files to metadata
###########################################################

mapping <- data.frame(snv = colnames(mat))
mapping$Sample.Name <- sub("^([^-]+-[^-]+-[^-]+).*", "\\1", mapping$snv)
mapping$Sample.Name <- gsub("-", "\\.", mapping$Sample.Name)

meta$snv_label <- mapping$snv[match(meta$Sample.Name, mapping$Sample.Name)]

###########################################################
# Formating dataframe for plotting
###########################################################

# list mutations of interest
common_mut <- c("PIK3CA", "TP53", "CDH1", "ERBB2", "MUC16", "TG", "TTN",
                 "QSER1", "SPTA1", "PTPRB", "GATA3", "USH2A", "TMCC3", "NCOR1", "ZFHX4", 
                 "PTEN", "BRCA1", "BRCA2", "ATM", "CHEK2")
toPlot <- mat[rownames(mat) %in% common_mut,]

# save labels for plotting
toPlot <- reshape2::melt(toPlot)
toPlot$Signature <- meta$Signature[match(toPlot$Var2, meta$snv_label)]
toPlot$Subtype <- meta$Subtype[match(toPlot$Var2, meta$snv_label)]
toPlot$rank <- meta$rank[match(toPlot$Var2, meta$snv_label)]

# format for plotting
toPlot <- toPlot[order(toPlot$rank),]
toPlot$rank <- factor(toPlot$rank, levels = c(75:1))
toPlot$value <- ifelse(toPlot$value >= 1, 1, 0)
toPlot$value <- factor(toPlot$value, levels = c(1, 0))

# format annotation bars
toPlot_map <- toPlot[,colnames(toPlot) %in% c("Var2", "Signature", "Subtype", "rank")]
#toPlot_map <- reshape2::melt(toPlot_map, id = c("Var2", "rank"))
toPlot_map <- toPlot_map[order(toPlot_map$rank),]
toPlot_map$rank <- factor(toPlot_map$rank, levels = c(75:1))

###########################################################
# Plot detection of common BCa mutations
###########################################################

# set up palette for plotting
pal <- c("Signature1" = "#046C9A", "Signature2" = "#BBADB9", "Signature3" = "#7294D4", 
         "Signature4" = "#E8E1D9", "Signature5" = "#AFC5D8", "Signature6" = "#DF9C93",
         "Signature7" = "#4E4C67", "Signature8" = "#B4869F", "Signature9" = "#A6B1E1",
         "Signature10" = "#985F6F")
#subtype_pal <- c("Basal" = "#394032","Her2" = "#A6A57A", "LumA" = "#8C5E58", "LumB" = "#5A352A", "Normal" = "#8F8073", "Not Available" = "#eFeBF7")
subtype_pal <- c("Basal" = "#AF4C5B","Her2" = "#EED4D3", "LumA" = "#B3B4D0", "LumB" = "#363E62", "Normal" = "#6365AF", "Not Available" = "#eFeBF7")


# plot heatmap
p1 <- ggplot(toPlot) + geom_tile(aes(x = Var1, y = rank, fill = value), color = "gray") +
  scale_fill_manual("Mutation Status", values = c("#533B4D", "white"), labels = c("Mutated", "Not Mutated")) +
  geom_hline(yintercept = c(16.5, 24.5, 40.5, 49.5, 58.5), linetype = "dotted", color = "gray") +
  theme_void() + 
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
  labs(x = "Genes", y = "Tumour Sample\n")

# signature annotation bar
p2 <- ggplot(toPlot_map) + geom_tile(aes(x = 1, y = rank, fill = Signature)) +
  theme_void() +
  scale_fill_manual(values = pal) +
  theme(axis.title.x = element_text(size=12, angle = 90, vjust = 0.5)) + 
  labs(fill = "Signature          ", x = "              ")

# subtype annotation bar
p3 <- ggplot(toPlot_map) + geom_tile(aes(x = 1, y = rank, fill = Subtype)) +
  theme_void() +
  scale_fill_manual(values = subtype_pal) +
  theme(axis.title.x = element_text(size=12, angle = 90, vjust = 0.5)) + 
  labs(fill = "Subtype      ", x = "              ")


# extract legends
l1 <- as_ggplot(get_legend(p1))
l2 <- as_ggplot(get_legend(p2))
l3 <- as_ggplot(get_legend(p3))
p1 <- p1+theme(legend.position = "none")
p2 <- p2+theme(legend.position = "none")
p3 <- p3+theme(legend.position = "none")

png("MolecularSigAnalysis/results/figures/heatmap.png", width = 8, height = 6, res = 600, units = "in")
grid.arrange(p1, p2, p3, l1, l2, l3, ncol = 19, nrow = 6,
             layout_matrix = rbind(c(1,1,1,1,1,1,1,1,1,1,1,2,3,4,4,4,4), 
                                   c(1,1,1,1,1,1,1,1,1,1,1,2,3,5,5,5,5), 
                                   c(1,1,1,1,1,1,1,1,1,1,1,2,3,5,5,5,5), 
                                   c(1,1,1,1,1,1,1,1,1,1,1,2,3,6,6,6,6),
                                   c(1,1,1,1,1,1,1,1,1,1,1,2,3,6,6,6,6),
                                   c(1,1,1,1,1,1,1,1,1,1,1,2,3,NA,NA,NA,NA)))
dev.off()