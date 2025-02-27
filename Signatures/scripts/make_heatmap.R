setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(data.table)
    library(ggplot2)
    library(reshape2)
    library(dplyr)
    library(RColorBrewer)
    library(ggh4x)
    library(ggpubr)
    library(grid)
    library(gridExtra)
})


###########################################################
# Load in data
###########################################################

# read in meta data
meta <- read.csv("Signatures/data/TCGA_subtype_label.csv")

# load in matrix file from NMF
mat <- read.table("Signatures/results/data/ATAC_heatmap_rank6.png.order.matrix", header = T)


###########################################################
# Format signature data for plotting
###########################################################

# get signatures
mat$Signature <- paste0("Signature", 1:6)

# format for plotting
mat <- melt(mat)

# save molecular subtype
meta$File.Name <- gsub("-", "\\.", meta$File.Name)
subtype <- c()
for (sample in gsub("X", "", unique(mat$variable))) {subtype <- c(subtype, meta[meta$File.Name == sample,]$bca_subtype)}
mat$subtype <- rep(subtype, each = 6)

# save signature assignment
mat$signature_assign <- ""
for (sample in mat$variable) {
    tmp <- mat[mat$variable == sample,]
    mat[mat$variable == sample,]$signature_assign <- as.character(tmp[which.max(tmp$value),]$Signature)
}

# reorder signatures
mat$Signature <- paste0(mat$Signature, " ")
mat$Signature <- factor(mat$Signature, levels = paste0("Signature",6:1," "))


###########################################################
# Set up palettes for heatmaps
###########################################################

# signature palette
pal <- c("Signature1" = "#046C9A", "Signature2" = "#BBADB9", "Signature3" = "#7294D4", 
        "Signature4" = "#E8E1D9", "Signature5" = "#AFC5D8", "Signature6" = "#DF9C93",
        "Signature7" = "#4E4C67", "Signature8" = "#B4869F", "Signature9" = "#A6B1E1",
        "Signature10" = "#985F6F")
# subtype palette
subtype_pal <- c("Basal" = "#AF4C5B","Her2" = "#EED4D3", "LumA" = "#B3B4D0", "LumB" = "#363E62", "Normal" = "#6365AF", "Not Available" = "#eFeBF7")


###########################################################
# Plot ATAC-Signature heatmaps
###########################################################

p1 <- ggplot(mat, aes(x = variable, y = 1, fill = subtype)) + geom_tile(color = NA) +
    theme_void() + scale_fill_manual("BCa\nSubtype", values = subtype_pal,
                    labels = c("Basal", "Her2", "LuminalA", "LuminalB", "Normal", "NA")) +
    theme(axis.title.y = element_text(size=12)) + labs(y = "                 ")

p2 <- ggplot(mat, aes(x = variable, y = 1, fill = signature_assign)) + geom_tile(color = NA) +
    theme_void() + scale_fill_manual("Assigned\nSignature", values = pal) +
    theme(axis.title.y = element_text(size=12)) + labs(y = "                 ")

p3 <- ggplot(mat, aes(x = variable, y = Signature, fill = value)) + geom_tile(color = NA) +
    scale_fill_gradientn("Signature    \nScore", colours = brewer.pal(9, "Blues")) + theme_void() +
    theme(axis.text.y = element_text(size=11), axis.title.x = element_text(size=12)) + labs(x = "Tumour Sample")

# extract legends
l1 <- as_ggplot(get_legend(p1))
l2 <- as_ggplot(get_legend(p2))
l3 <- as_ggplot(get_legend(p3))
p1 <- p1+theme(legend.position = "none")
p2 <- p2+theme(legend.position = "none")
p3 <- p3+theme(legend.position = "none")

png("Signatures/results/figures/ATAC_rank6_heatmap.png", width = 11, height = 4, res = 600, units = "in")
grid.arrange(p1, p2, p3, l1, l2, l3, ncol = 9, nrow = 12,
    layout_matrix = rbind(c(1,1,1,1,1,1,1,4,5), c(2,2,2,2,2,2,2,4,5),
                          c(3,3,3,3,3,3,3,4,5), c(3,3,3,3,3,3,3,4,5),
                          c(3,3,3,3,3,3,3,4,5), c(3,3,3,3,3,3,3,4,5),
                          c(3,3,3,3,3,3,3,6,NA),c(3,3,3,3,3,3,3,6,NA),
                          c(3,3,3,3,3,3,3,6,NA),c(3,3,3,3,3,3,3,6,NA),
                          c(3,3,3,3,3,3,3,6,NA),c(3,3,3,3,3,3,3,NA,NA)))
dev.off()


###########################################################
# Format RNA-Seq NMF matrices for heatmap
###########################################################

# function for processing nmf heatmap files for RNA-Seq
read_mat <- function(file, rank) {
        
    # load in matrix file from NMF
    mat <- read.table(file, header = T)

    # get signature
    mat$Signature <- paste0("Signature", 1:rank)

    # format for plotting
    mat <- melt(mat)

    # save molecular subtype
    subtype <- c()
    for (sample in mat$variable) {subtype <- c(subtype, meta[meta$sample_name == sample,]$bca_subtype[1])}
    mat$subtype <- subtype

    # save signature assignment
    mat$signature_assign <- ""
    for (sample in mat$variable) {
        tmp <- mat[mat$variable == sample,]
        mat[mat$variable == sample,]$signature_assign <- as.character(tmp[which.max(tmp$value),]$Signature)
    }

    # reorder signatures
    mat$Signature <- factor(mat$Signature, levels = paste0("Signature",rank:1))

    return(mat)
}

mat <- read_mat("Signatures/results/data/RNA_heatmap_rank6.png.order.matrix", 6)


###########################################################
# Plot RNA-Signature heatmap
###########################################################

p1 <- ggplot(mat, aes(x = variable, y = 1, fill = subtype)) + geom_tile(color = NA) +
    theme_void() + scale_fill_manual("BCa\nSubtype", values = subtype_pal,
                    labels = c("Basal", "Her2", "LuminalA", "LuminalB", "Normal", "NA")) +
    theme(axis.title.y = element_text(size=12)) + labs(y = "                 ")

p2 <- ggplot(mat, aes(x = variable, y = 1, fill = signature_assign)) + geom_tile(color = NA) +
    theme_void() + scale_fill_manual("Assigned\nSignature", values = pal) +
    theme(axis.title.y = element_text(size=12)) + labs(y = "                 ")

p3 <- ggplot(mat, aes(x = variable, y = Signature, fill = value)) + geom_tile(color = NA) +
    scale_fill_gradientn("Signature    \nScore", colours = brewer.pal(9, "Blues")) + theme_void() +
    theme(axis.text.y = element_text(size=11), axis.title.x = element_text(size=12)) + labs(x = "Tumour Sample")

# extract legends
l1 <- as_ggplot(get_legend(p1))
l2 <- as_ggplot(get_legend(p2))
l3 <- as_ggplot(get_legend(p3))
p1 <- p1+theme(legend.position = "none")
p2 <- p2+theme(legend.position = "none")
p3 <- p3+theme(legend.position = "none")

png("Signatures/results/figures/RNA_rank6_heatmap.png", width = 11, height = 4, res = 600, units = "in")
grid.arrange(p1, p2, p3, l1, l2, l3, ncol = 9, nrow = 12,
    layout_matrix = rbind(c(1,1,1,1,1,1,1,4,5), c(2,2,2,2,2,2,2,4,5),
                          c(3,3,3,3,3,3,3,4,5), c(3,3,3,3,3,3,3,4,5),
                          c(3,3,3,3,3,3,3,4,5), c(3,3,3,3,3,3,3,4,5),
                          c(3,3,3,3,3,3,3,6,NA),c(3,3,3,3,3,3,3,6,NA),
                          c(3,3,3,3,3,3,3,6,NA),c(3,3,3,3,3,3,3,6,NA),
                          c(3,3,3,3,3,3,3,6,NA),c(3,3,3,3,3,3,3,NA,NA)))
dev.off()
