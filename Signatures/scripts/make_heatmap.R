setwd("C:/Users/julia/Documents/BCaATAC")

library(data.table)
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(ggh4x)
library(ggpubr)
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

setwd("C:/Users/julia/Documents/BCaATAC")


# set up palette for plotting
pal <- c("Signature1" = "#046C9A", "Signature2" = "#BBADB9", "Signature3" = "#7294D4", 
        "Signature4" = "#E8E1D9", "Signature5" = "#AFC5D8", "Signature6" = "#DF9C93",
        "Signature7" = "#4E4C67", "Signature8" = "#B4869F", "Signature9" = "#A6B1E1",
        "Signature10" = "#985F6F")
subtype_pal <- c("#394032", "#A6A57A", "#8C5E58", "#5A352A", "#8F8073", "#eFeBF7")

# load in matrix file from NMF
mat <- read.table("Signatures/results/ATAC_heatmap_rank6.png.order.matrix", header = T)

# get signature
mat$Signature <- paste0("Signature", 1:6)

# format for plotting
mat <- melt(mat)

# read in meta data
meta <- fread("Signatures/data/brca.assignment.share.txt")
meta <- meta[1:75]

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
mat$Signature <- factor(mat$Signature, levels = paste0("Signature",6:1))



### Plot Heatmap Rank 6 ###
p1 <- ggplot(mat, aes(x = variable, y = 1, fill = subtype)) + geom_tile(color = "black") +
    theme_void() + scale_fill_manual("BCa\nSubtype", values = subtype_pal,
                    labels = c("Basal", "Her2", "LuminalA", "LuminalB", "Normal", "NA")) +
    theme(axis.title.y = element_text(size=12)) + labs(y = "                ")

p2 <- ggplot(mat, aes(x = variable, y = 1, fill = signature_assign)) + geom_tile(color = "black") +
    theme_void() + scale_fill_manual("Assigned\nSignature", values = pal) +
    theme(axis.title.y = element_text(size=12)) + labs(y = "                ")

p3 <- ggplot(mat, aes(x = variable, y = Signature, fill = value)) + geom_tile(color = "black") +
    scale_fill_gradientn("Signature    \nScore", colours = brewer.pal(9, "Blues")) + theme_void() +
    theme(axis.text.y = element_text(size=11), axis.title.x = element_text(size=12)) + labs(x = "Tumour Sample")

# extract legends
l1 <- as_ggplot(get_legend(p1))
l2 <- as_ggplot(get_legend(p2))
l3 <- as_ggplot(get_legend(p3))
p1 <- p1+theme(legend.position = "none")
p2 <- p2+theme(legend.position = "none")
p3 <- p3+theme(legend.position = "none")

png("Signatures/results/rank6_heatmap.png", width = 10, height = 4, res = 600, units = "in")
grid.arrange(p1, p2, p3, l1, l2, l3, ncol = 9, nrow = 13,
    layout_matrix = rbind(c(1,1,1,1,1,1,1,4,5),
                        c(2,2,2,2,2,2,2,4,5),
                        c(3,3,3,3,3,3,3,4,5),
                        c(3,3,3,3,3,3,3,4,5),
                        c(3,3,3,3,3,3,3,4,5),
                        c(3,3,3,3,3,3,3,4,5),
                        c(3,3,3,3,3,3,3,4,5),
                        c(3,3,3,3,3,3,3,6,NA),
                        c(3,3,3,3,3,3,3,6,NA),
                        c(3,3,3,3,3,3,3,6,NA),
                        c(3,3,3,3,3,3,3,6,NA),
                        c(3,3,3,3,3,3,3,6,NA),
                        c(3,3,3,3,3,3,3,NA,NA)))
dev.off()