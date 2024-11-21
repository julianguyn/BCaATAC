# Script to investigate feature importance of ElasticNet and LASSO

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(GenomicRanges)
    library(RColorBrewer)
    library(tidyverse)
    library(reshape2)
})

# load in annotation 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# define function for peak annotation
peak_anno <- function(features) {
    # extract feature information
    features <- features %>%
        separate(Peak, into = c("chr", "start", "end"), sep = ":", convert = TRUE)

    # create GRanges object
    gr <- GRanges(seqnames = features$chr, ranges = IRanges(features$start, features$end))

    # annotate peaks
    peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

    return(peakAnno)
}

### ===== Feature Importance of Elastic Net ===== ###

# load in feature importance
en_feat <- read.csv("SupervisedSignatures/results/data/feature_importance_elasticnet.csv")

# get number of excluded features
nrow(en_feat[en_feat$Weight == 0,])
# 981 out of 99

# subset for important features (left with 18 features)
en_feat <- en_feat[en_feat$Weight != 0,]

# formatting for plotting
en_feat$Peak <- factor(en_feat$Peak, levels = en_feat$Peak)

# plot feature importance
png("SupervisedSignatures/results/figures/elasticnet_features.png", width=200, height=100, units='mm', res = 600, pointsize=80)
ggplot(en_feat, aes(x = Peak, y = abs(Weight), fill = ifelse(Weight > 0, "Positive", "Negative"))) + 
    geom_bar(stat="identity") + theme_classic() +
    scale_fill_manual(values = c("#7294D4", "#BC4749"), limits = c("Positive", "Negative")) +
    labs(fill = "Association", y = "Coefficient Magnitude") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.key.size = unit(0.7, 'cm')) 
dev.off()


# plot peak distribution
png("SupervisedSignatures/results/figures/elasticnet_peakAnno.png", width = 6, height = 5, res = 600, units = "in")
plotAnnoPie(peak_anno(en_feat))
dev.off()



### ===== Feature Importance of LASSO ===== ###

# load in feature importance
ls_feat <- read.csv("SupervisedSignatures/results/data/feature_importance_lasso.csv")

# get number of excluded features
nrow(ls_feat[ls_feat$Weight == 0,])
# 989 out of 99

# subset for important features (left with 10 features)
ls_feat <- ls_feat[ls_feat$Weight != 0,]

# formatting for plotting
ls_feat$Peak <- factor(ls_feat$Peak, levels = ls_feat$Peak)

# plot feature importance
png("SupervisedSignatures/results/figures/lasso_features.png", width=200, height=100, units='mm', res = 600, pointsize=80)
ggplot(ls_feat, aes(x = Peak, y = abs(Weight), fill = ifelse(Weight > 0, "Positive", "Negative"))) + 
    geom_bar(stat="identity") + theme_classic() +
    scale_fill_manual(values = c("#7294D4", "#BC4749"), limits = c("Positive", "Negative")) +
    labs(fill = "Association", y = "Coefficient Magnitude") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.key.size = unit(0.7, 'cm')) 
dev.off()


# plot peak distribution
png("SupervisedSignatures/results/figures/lasso_peakAnno.png", width = 6, height = 5, res = 600, units = "in")
plotAnnoPie(peak_anno(ls_feat))
dev.off()