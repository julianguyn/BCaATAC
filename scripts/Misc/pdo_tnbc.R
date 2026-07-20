# for parisa (lupien lab)
# get the binary matrix

library(data.table)


bin <- readRDS("data/rawdata/misc/pdo_tnbc.consensus.Binarymat.rds")
colnames(bin)[1:3] <- c("V1", "V2", "V3")
write.table(bin, file = "data/rawdata/misc/pdo_tnbc.matrix", sep = "\t", quote = FALSE, row.names = FALSE)


df <- fread('data/rawdata/tcga/BCa_binary.2.matrix')

df <- fread("data/rawdata/misc/pdo_tnbc.matrix")
