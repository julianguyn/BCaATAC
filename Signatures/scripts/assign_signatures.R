setwd("C:/Users/julia/Documents/BCaATAC")
mat <- read.table("Signatures/bca_sign.Zscore.txt")

signatures <- c() 
scores <- c()
for (s in 1:ncol(mat)) {
    # check for ties
    max_score <- max(mat[,s])
    if (max_score %in% mat[,s][duplicated(mat[,s])]) { print("tie") }

    # save signature with the maximum chromvar score and the score
    signatures <- c(signatures, rownames(mat)[which.max(mat[,s])])
    scores <- c(scores, max_score)
}

df <- data.frame(sample = colnames(mat), signature = signatures, score = scores)

order <- read.table("Signatures/chromvar_order.txt")
df <- df[match(order$V1, df$sample),]

write.table(df, file = "Signatures/signature_scores.tsv", quote = F, sep = "\t", col.names = T, row.names = F)