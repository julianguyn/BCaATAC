### Script to compare overlap in TCGA samples between ATAC and RNA signatures


setwd("C:/Users/julia/Documents/BCaATAC")

# set up palette for plotting
pal <- c("Signature1" = "#046C9A", "Signature2" = "#BBADB9", "Signature3" = "#7294D4", 
        "Signature4" = "#E8E1D9", "Signature5" = "#AFC5D8", "Signature6" = "#DF9C93",
        "Signature7" = "#4E4C67", "Signature8" = "#B4869F", "Signature9" = "#A6B1E1",
        "Signature10" = "#985F6F")
subtype_pal <- c("#394032", "#A6A57A", "#8C5E58", "#5A352A", "#8F8073", "#eFeBF7")

# read in meta data
meta <- read.csv("Signatures/data/TCGA_subtype_label.csv")
meta$File.Name <- gsub("-", "\\.", meta$File.Name)

# function to read in binary matrix files
read_mat <- function(file, rank, atac=FALSE) {
        
    # load in matrix file from NMF
    mat <- read.table(file, header = T)

    # get signature
    mat$Signature <- paste0("Signature", 1:rank)

    # format for plotting
    mat <- melt(mat)

    # save signature assignment
    mat$signature_assign <- ""
    for (sample in mat$variable) {
        tmp <- mat[mat$variable == sample,]
        mat[mat$variable == sample,]$signature_assign <- as.character(tmp[which.max(tmp$value),]$Signature)
    }

    if (atac == TRUE) {

        # reorder signatures
        mat$Signature <- factor(mat$Signature, levels = paste0("Signature",rank:1))

        # match names to RNA-Seq
        mat$sample_name <- ""
        for (i in 1:nrow(mat)) {mat$sample_name[i] <- meta[which(meta$File.Name == gsub("X", "", mat$variable[i])),]$sample_name}

        # keep unique sample - signature assignments
        mat <- unique(mat[,c(4:5)])
    } else {
        # keep unique sample - signature assignments
        mat <- unique(mat[,c(4,2)])
    }

    return(mat)
}


# load in ATAC matrix file from NMF
atac_mat <- read_mat("Signatures/results/data/ATAC_heatmap_rank6.png.order.matrix", 6, TRUE)

# function to process RNA matrix file and get number of overlapping samples
get_overlap <- function(filename, rank) {

    # read in RNA matrix file from NMF
    rna_mat <- read_mat(filename, rank)

    # get all signature combinations
    combinations <- expand.grid(unique(atac_mat$signature_assign), unique(rna_mat$signature_assign))
    colnames(combinations) <- c("ATACSeq", "RNASeq")
    combinations$count <- 0

    # get number of overlapping tumour samples
    for (i in 1:nrow(combinations)) {

        # subset signatures from atac and rna
        atac_s <- atac_mat[atac_mat$signature_assign == combinations$ATACSeq[i],]
        rna_s <- rna_mat[rna_mat$signature_assign == combinations$RNASeq[i],]

        # get number of intersecting samples
        n_common <- length(intersect(atac_s$sample_name, rna_s$variable))

        # populate combinations with number
        combinations$count[i] <- n_common
    }

    # reorder RNA-Seq signatures
    combinations$RNASeq <- factor(combinations$RNASeq, levels = paste0("Signature",rank:1))

    p <- ggplot(combinations, aes(x = ATACSeq, y = RNASeq, fill = count)) + geom_tile(color = "black") +
            scale_fill_gradientn("No. Overlapping\nSamples", colours = brewer.pal(9, "Blues")) + theme_void() +
            theme(axis.text.x = element_text(size=11, angle = 90, vjust = 0.5), axis.title.x = element_text(size=12),
                axis.text.y = element_text(size=11), axis.title.y = element_text(size=12, angle = 90, vjust = 0.5)) + 
            labs(x = "\nATAC-Seq Signatures", y = "RNA-Seq Signatures\n")

    return(p)
}

png("Signatures/results/figures/overlap_rank4.png", width = 4, height = 6, res = 600, units = "in")
get_overlap("Signatures/results/data/RNA_heatmap_rank4.png.order.matrix", 4)
dev.off()

png("Signatures/results/figures/overlap_rank5.png", width = 4, height = 6, res = 600, units = "in")
get_overlap("Signatures/results/data/RNA_heatmap_rank5.png.order.matrix", 5)
dev.off()

png("Signatures/results/figures/overlap_rank6.png", width = 4, height = 6, res = 600, units = "in")
get_overlap("Signatures/results/data/RNA_heatmap_rank6.png.order.matrix", 6)
dev.off()

png("Signatures/results/figures/overlap_rank7.png", width = 4, height = 6, res = 600, units = "in")
get_overlap("Signatures/results/data/RNA_heatmap_rank7.png.order.matrix", 7)
dev.off()

png("Signatures/results/figures/overlap_rank8.png", width = 4, height = 6, res = 600, units = "in")
get_overlap("Signatures/results/data/RNA_heatmap_rank8.png.order.matrix", 8)
dev.off()

png("Signatures/results/figures/overlap_rank9.png", width = 4, height = 6, res = 600, units = "in")
get_overlap("Signatures/results/data/RNA_heatmap_rank9.png.order.matrix", 9)
dev.off()