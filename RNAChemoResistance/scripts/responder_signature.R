setwd("C:/Users/julia/Documents/BCaATAC")

library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(ggh4x)
library(ggpubr)suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

# set up palette for plotting
pal <- c("#899DA4", "#BC4749")

# load in signature scores
rank5 <- read.table("RNAChemoResistance/results/data/heatmap_rank5.png.order.matrix", header = T)
rank6 <- read.table("RNAChemoResistance/results/data/heatmap_rank6.png.order.matrix", header = T)


# extract signatures
rank5$Signature <- paste0("Signature", 1:5)
rank6$Signature <- paste0("Signature", 1:6)
rank5 <- melt(rank5)
rank6 <- melt(rank6)

# load in meta data
meta <- read.csv("RNAChemoResistance/data/metadata.csv")
meta$label <- paste0(meta$Sampleid, meta$Rep)

# get responder status
rank5$responder <- meta[match(rank5$variable, meta$label),]$Responder.status
rank6$responder <- meta[match(rank6$variable, meta$label),]$Responder.status
rank5$responder <- factor(rank5$responder, levels = c(TRUE, FALSE))
rank6$responder <- factor(rank6$responder, levels = c(TRUE, FALSE))

# modify signature levels
rank5$Signature <- factor(rank5$Signature, levels = paste0("Signature", 5:1))
rank6$Signature <- factor(rank6$Signature, levels = paste0("Signature", 6:1))


### Plot Heatmap Rank 5 ###
p1 <- ggplot(rank5, aes(x = variable, y = 1, fill = responder)) + geom_tile(color = "black") +
    theme_void() + scale_fill_manual("Responder\nScore", values = pal) +
    theme(axis.title.y = element_text(size=12)) + labs(y = "                ")

p2 <- ggplot(rank5, aes(x = variable, y = Signature, fill = value)) + geom_tile(color = "black") +
    scale_fill_gradientn("Signature \nScore", colours = brewer.pal(9, "Blues")) + theme_void() +
    theme(axis.text.y = element_text(size=11))

# extract legends
l1 <- as_ggplot(get_legend(p1))
l2 <- as_ggplot(get_legend(p2))
p1 <- p1+theme(legend.position = "none")
p2 <- p2+theme(legend.position = "none")

png("RNAChemoResistance/results/figures/rank5_heatmap.png", width = 8, height = 3, res = 600, units = "in")
grid.arrange(p1, p2, l1, l2, ncol = 8, nrow = 8,
    layout_matrix = rbind(c(1,1,1,1,1,1,1,3),
                        c(2,2,2,2,2,2,2,3),
                        c(2,2,2,2,2,2,2,3),
                        c(2,2,2,2,2,2,2,3),
                        c(2,2,2,2,2,2,2,4),
                        c(2,2,2,2,2,2,2,4),
                        c(2,2,2,2,2,2,2,4),
                        c(2,2,2,2,2,2,2,NA)))
dev.off()


### Plot Heatmap Rank 6 ###
p1 <- ggplot(rank6, aes(x = variable, y = 1, fill = responder)) + geom_tile(color = "black") +
    theme_void() + scale_fill_manual("Responder\nScore", values = pal) +
    theme(axis.title.y = element_text(size=12)) + labs(y = "                ")

p2 <- ggplot(rank6, aes(x = variable, y = Signature, fill = value)) + geom_tile(color = "black") +
    scale_fill_gradientn("Signature \nScore", colours = brewer.pal(9, "Blues")) + theme_void() +
    theme(axis.text.y = element_text(size=11))

# extract legends
l1 <- as_ggplot(get_legend(p1))
l2 <- as_ggplot(get_legend(p2))
p1 <- p1+theme(legend.position = "none")
p2 <- p2+theme(legend.position = "none")

png("RNAChemoResistance/results/figures/rank6_heatmap.png", width = 8, height = 3, res = 600, units = "in")
grid.arrange(p1, p2, l1, l2, ncol = 8, nrow = 8,
    layout_matrix = rbind(c(1,1,1,1,1,1,1,3),
                        c(2,2,2,2,2,2,2,3),
                        c(2,2,2,2,2,2,2,3),
                        c(2,2,2,2,2,2,2,3),
                        c(2,2,2,2,2,2,2,4),
                        c(2,2,2,2,2,2,2,4),
                        c(2,2,2,2,2,2,2,4),
                        c(2,2,2,2,2,2,2,NA)))
dev.off()




### Plot Proportions ###

rank5 <- as.data.frame(rank5 %>% group_by(variable) %>% filter(value == max(value)) %>% ungroup())
rank6 <- as.data.frame(rank6 %>% group_by(variable) %>% filter(value == max(value)) %>% ungroup())


# plot distribution of responders by signature
png("RNAChemoResistance/results/figures/rank5_proportion.png", width = 6, height = 4, res = 600, units = "in")
ggplot(rank5, aes(x = Signature, fill = factor(responder))) + geom_bar(position = "fill") + 
    theme_classic() + labs(y = "Proportion") + scale_fill_manual("Responder Status", values = pal)
dev.off()
png("RNAChemoResistance/results/figures/rank6_proportion.png", width = 6, height = 4, res = 600, units = "in")
ggplot(rank6, aes(x = Signature, fill = factor(responder))) + geom_bar(position = "fill") + 
    theme_classic() + labs(y = "Proportion") + scale_fill_manual("Responder Status", values = pal)
dev.off()
