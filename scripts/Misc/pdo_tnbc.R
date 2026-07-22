# for parisa (lupien lab)
# get the binary matrix

library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)


bin <- readRDS("data/rawdata/misc/pdo_tnbc.consensus.Binarymat.rds")
colnames(bin)[1:3] <- c("V1", "V2", "V3")
write.table(bin, file = "data/rawdata/misc/pdo_tnbc.matrix", sep = "\t", quote = FALSE, row.names = FALSE)


df <- fread('data/rawdata/tcga/BCa_binary.2.matrix')

df <- fread("data/rawdata/misc/pdo_tnbc.matrix")


# make heatmap

df <- read.table("data/results/data/Misc/PDO_TNBC/rank5.png.order.matrix", header = T)
rownames(df) <- paste0("Signature", 1:nrow(df))

# make colour palette
cols <- brewer.pal(9, "Blues")
col_fun <- colorRamp2(
    seq(min(df, na.rm = TRUE),
        max(df, na.rm = TRUE),
        length.out = 9),
    cols
)

# heatmap
ht1 <- Heatmap(
    df,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "NMF\nWeight",
    col = col_fun,
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 8)
)

filename <- "data/results/figures/pdo_tnbc_rank5_heatmap.png"
png(filename, width = 5, height = 4, res = 600, units = "in")
ht1
dev.off()
