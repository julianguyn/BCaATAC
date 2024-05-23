# Script to look at the distribution of peak calls for all used samples

setwd("C:/Users/julia/Documents/BCaATAC")

suppressMessages(library(ComplexHeatmap))

# set up palette for plotting
pal <- c("Signature1" = "#046C9A", "Signature2" = "#BBADB9", "Signature3" = "#7294D4", 
        "Signature4" = "#E8E1D9", "Signature5" = "#AFC5D8", "Signature6" = "#DF9C93",
        "Signature7" = "#4E4C67", "Signature8" = "#B4869F", "Signature9" = "#A6B1E1",
        "Signature10" = "#985F6F")
subtype_pal <- c("#394032", "#A6A57A", "#8C5E58", "#5A352A", "#8F8073", "#eFeBF7")

### CHECK PEAK DISTRIBUTION OF TUMOURS

# load in TCGA binary matrix
TCGA <- read.table("Signatures/data/BCa_binary.2.matrix", header = T)
TCGA$peaks <- paste(TCGA$V1, TCGA$V2, TCGA$V3, sep = "-")
TCGA <- TCGA[,-c(1:3)]

# read in meta data
meta <- read.csv("Signatures/data/TCGA_subtype_label.csv")

# load in matrix file from NMF
mat <- read.table("Signatures/results/data/ATAC_heatmap_rank6.png.order.matrix", header = T)

# get signature
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

# save labels
labels <- unique(mat[,c(2,4,5)])
labels <- labels[match(colnames(TCGA)[1:75], labels$variable),]


# format dataframe for plotting
#npeaks <- nrow(TCGA)
#TCGA <- melt(TCGA)
#TCGA$subtype <- rep(labels$subtype, each = npeaks)
#TCGA$signature <- rep(labels$signature_assign, each = npeaks)

library(circlize)
sig_score = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
sig_score(seq(-3, 3))

m = t(TCGA[,c(1:75)])
pdf(file = "DataExploration/results/figures/heatmap.pdf")
Heatmap(m, name = "Signature Score", col = sig_score)
dev.off()





ht_list = Heatmap(m, name = "main_matrix")




ha = HeatmapAnnotation(summary = anno_summary(height = unit(3, "cm")))
ht_list = ht_list + Heatmap(labels$subtype, name = "subtype", top_annotation = ha, width = unit(1, "cm"))

ha = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:3), 
    height = unit(3, "cm")))
ht_list = ht_list + Heatmap(labels$signature_assign, name = "signature", top_annotation = ha, width = unit(1, "cm"))

split = labels$signature_assign
lgd_boxplot = Legend(labels = names(pal), title = "Signature",
    legend_gp = gpar(fill = pal))
draw(ht_list, row_split = split, ht_gap = unit(5, "mm"), 
    heatmap_legend_list = list(lgd_boxplot))




# original code:
m = matrix(rnorm(50*10), nrow = 50)
ht_list = Heatmap(m, name = "main_matrix")

ha = HeatmapAnnotation(summary = anno_summary(height = unit(3, "cm")))
v = sample(letters[1:2], 50, replace = TRUE)
ht_list = ht_list + Heatmap(v, name = "mat1", top_annotation = ha, width = unit(1, "cm"))

ha = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:3), 
    height = unit(3, "cm")))
v = rnorm(50)
ht_list = ht_list + Heatmap(v, name = "mat2", top_annotation = ha, width = unit(1, "cm"))

split = sample(letters[1:2], 50, replace = TRUE)
lgd_boxplot = Legend(labels = c("group a", "group b"), title = "group",
    legend_gp = gpar(fill = c("red", "blue")))
draw(ht_list, row_split = split, ht_gap = unit(5, "mm"), 
    heatmap_legend_list = list(lgd_boxplot))


# plot distribution of peaks
p1 <- ggplot(TCGA, aes(x = 1, y = variable, fill = subtype)) + geom_tile(color = "black") +
    theme_void() + scale_fill_manual("BCa\nSubtype", values = subtype_pal,
                    labels = c("Basal", "Her2", "LuminalA", "LuminalB", "Normal", "NA")) +
    theme(axis.title.y = element_text(size=12)) + labs(y = "                ")

p2 <- ggplot(TCGA, aes(x = 1, y = variable, fill = signature_assign)) + geom_tile(color = "black") +
    theme_void() + scale_fill_manual("Assigned\nSignature", values = pal) +
    theme(axis.title.y = element_text(size=12)) + labs(y = "                ")

p3 <- ggplot(TCGA, aes(x = peaks, y = variable, fill = value)) + geom_tile(color = "black") +
    scale_fill_gradientn("Signature    \nScore", colours = brewer.pal(9, "Blues")) + theme_void() +
    theme(axis.text.y = element_text(size=11), axis.title.x = element_text(size=12)) + labs(x = "Tumour Sample")

# extract legends
l1 <- as_ggplot(get_legend(p1))
l2 <- as_ggplot(get_legend(p2))
l3 <- as_ggplot(get_legend(p3))
p1 <- p1+theme(legend.position = "none")
p2 <- p2+theme(legend.position = "none")
p3 <- p3+theme(legend.position = "none")

png("DataExploration/results/figures/tumour_peaks.png", width = 10, height = 6, res = 600, units = "in")
grid.arrange(p1, p2, p3, l1, l2, l3, ncol = 9, nrow = 13,
    layout_matrix = rbind(c(1,1,1,1,1,1,1,4,5),
                          c(1,2,2,2,2,2,2,4,5),
                          c(1,3,3,3,3,3,3,4,5),
                          c(1,3,3,3,3,3,3,4,5),
                          c(1,3,3,3,3,3,3,4,5),
                          c(1,3,3,3,3,3,3,4,5),
                          c(1,3,3,3,3,3,3,4,5),
                          c(1,3,3,3,3,3,3,3,6,NA),
                          c(1,3,3,3,3,3,3,6,NA),
                          c(1,3,3,3,3,3,3,6,NA),
                          c(1,3,3,3,3,3,3,6,NA),
                          c(1,3,3,3,3,3,3,6,NA),
                          c(1,3,3,3,3,3,3,NA,NA)))
dev.off()


### CHECK PEAK DISTRIBUTION OF CELL LINES

# load in cell line binary matrix
cells <- read.table("Signatures/data/bcacell_lines.tsv", header = T)

