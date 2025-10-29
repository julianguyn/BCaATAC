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
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(ComplexHeatmap)
    library(GenomicRanges)
    library(rGREAT)
    library(RColorBrewer)
})

source("utils/plots/signatures.R")
source("utils/get_data.R")
source("utils/palettes.R")

set.seed(123)


#' For now, using the 20k sites
#' May be worth assessing with 10k?
# analysis = c("20k", "all") 
analysis <- "20k"

###########################################################
# Load in data
###########################################################

# get siganture scores from NMF
mat <- get_arche_tcga()

###########################################################
# Plot ATAC-Signature heatmap
###########################################################

plot_ARCHE_heatmap(mat)

###########################################################
# Load in BEDs
###########################################################

bg <- fread(paste0("data/procdata/ARCHEs/beds/Background.bed"))

sig1 <- fread(paste0("data/procdata/ARCHEs/beds/ARCHE1_", analysis, ".bed"))
sig2 <- fread(paste0("data/procdata/ARCHEs/beds/ARCHE2_", analysis, ".bed"))
sig3 <- fread(paste0("data/procdata/ARCHEs/beds/ARCHE3_", analysis, ".bed"))
sig4 <- fread(paste0("data/procdata/ARCHEs/beds/ARCHE4_", analysis, ".bed"))
sig5 <- fread(paste0("data/procdata/ARCHEs/beds/ARCHE5_", analysis, ".bed"))
sig6 <- fread(paste0("data/procdata/ARCHEs/beds/ARCHE6_", analysis, ".bed"))

###########################################################
# Compute number and size of peak regions
###########################################################

# compute total open peak region
sig1$diff <- sig1$chromEnd - sig1$chromStart
sig2$diff <- sig2$chromEnd - sig2$chromStart
sig3$diff <- sig3$chromEnd - sig3$chromStart
sig4$diff <- sig4$chromEnd - sig4$chromStart
sig5$diff <- sig5$chromEnd - sig5$chromStart
sig6$diff <- sig6$chromEnd - sig6$chromStart

# get peak information
num_windows <- c(nrow(sig1), nrow(sig2), nrow(sig3), nrow(sig4), nrow(sig5), nrow(sig6))
sum_peaks <- c(sum(sig1$diff), sum(sig2$diff), sum(sig3$diff), sum(sig4$diff), sum(sig5$diff), sum(sig6$diff))

# plot peak info
df <- data.frame(ARCHE = paste0("ARCHE", 1:6),
                num_windows = num_windows, sum_peaks = sum_peaks)
plot_ARCHE_peakInfo(df, analysis)

###########################################################
# Compute number of overlapping regions
###########################################################

# create GRanges
grb <- GRanges(seqnames = paste0("chr", bg$chrom), ranges = IRanges(bg$chromStart, bg$chromEnd))

gr1 <- GRanges(seqnames = paste0("chr", sig1$chrom), ranges = IRanges(sig1$chromStart, sig1$chromEnd))
gr2 <- GRanges(seqnames = paste0("chr", sig2$chrom), ranges = IRanges(sig2$chromStart, sig2$chromEnd))
gr3 <- GRanges(seqnames = paste0("chr", sig3$chrom), ranges = IRanges(sig3$chromStart, sig3$chromEnd))
gr4 <- GRanges(seqnames = paste0("chr", sig4$chrom), ranges = IRanges(sig4$chromStart, sig4$chromEnd))
gr5 <- GRanges(seqnames = paste0("chr", sig5$chrom), ranges = IRanges(sig5$chromStart, sig5$chromEnd))
gr6 <- GRanges(seqnames = paste0("chr", sig6$chrom), ranges = IRanges(sig6$chromStart, sig6$chromEnd))

# create peak list for Upset plot
peak_list <- list(
    ARCHE1 = gr1, 
    ARCHE2 = gr2, 
    ARCHE3 = gr3, 
    ARCHE4 = gr4, 
    ARCHE5 = gr5, 
    ARCHE6 = gr6
)

# plot UPSET plot of overlapping peaks
m = make_comb_mat(peak_list)
plot_ATAC_Upset(m, analysis)

###########################################################
# Annotate ARCHE peak sets
###########################################################

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# helper function
annotateARCHE <- function(gr, arche) {
    anno <- annotatePeak(gr, tssRegion=c(-250, 250), TxDb=txdb, annoDb="org.Hs.eg.db")@annoStat
    anno$ARCHE <- arche
    return(anno)
}

# annotate peaks
anno1 <- annotateARCHE(gr1, "ARCHE1")
anno2 <- annotateARCHE(gr2, "ARCHE2")
anno3 <- annotateARCHE(gr3, "ARCHE3")
anno4 <- annotateARCHE(gr4, "ARCHE4")
anno5 <- annotateARCHE(gr5, "ARCHE5")
anno6 <- annotateARCHE(gr6, "ARCHE6")
annob <- annotateARCHE(grb, "Background")

# plot peakAnno results
toPlot <- rbind(anno1, anno2, anno3, anno4, anno5, anno6, annob)
plot_annotatePeak(toPlot, analysis)

###########################################################
# GREAT analysis
###########################################################

# helper function
runGREAT <- function(gr, arche, analysis) {
    
    job = submitGreatJob(gr, bg, species = "hg38", genome = "hg38", help = FALSE)
    tbl = getEnrichmentTables(job)

    # save results
    mf <- as.data.frame(tbl[1])
    bp <- as.data.frame(tbl[2])
    cc <- as.data.frame(tbl[3])

    # save column names
    cols <- gsub("GO.Molecular.Function.", "", colnames(mf))
    colnames(mf) <- colnames(bp) <- colnames(cc) <- cols

    mf$Label <- "Molecular Feature"
    bp$Label <- "Biological Process"
    cc$Label <- "Cellular Component"

    # combine results and filter
    res <- rbind(bp, mf, cc)
    res <- res[res$Hyper_Adjp_BH < 0.05,]
    res <- res[order(res$Hyper_Adjp_BH),]

    write.table(res, file = paste0("data/results/data/1-Signatures/GREAT/", arche, "_", analysis, ".tsv"), quote = F, sep = "\t", col.names = T, row.names = F)
    return(res)
}


# run GREAT
great1 <- runGREAT(gr1, "ARCHE1", analysis)
great2 <- runGREAT(gr2, "ARCHE2", analysis)
great3 <- runGREAT(gr3, "ARCHE3", analysis)
great4 <- runGREAT(gr4, "ARCHE4", analysis)
great5 <- runGREAT(gr5, "ARCHE5", analysis)
great6 <- runGREAT(gr6, "ARCHE6", analysis)

# great1 <- fread("data/results/data/1-Signatures/GREAT/ARCHE1_20k.tsv", data.table = FALSE)

###########################################################
# Plot GREAT enrichment
###########################################################

n <- 20

plot_GREAT(great1, n, "ARCHE1")
plot_GREAT(great2, n, "ARCHE2")
plot_GREAT(great3, n, "ARCHE3")
plot_GREAT(great4, n, "ARCHE4")
plot_GREAT(great5, n, "ARCHE5")
plot_GREAT(great6, n, "ARCHE6")

###########################################################
# Annotated GREAT enrichment plot
###########################################################

anno <- read.csv("data/procdata/ARCHEs/GREAT_20_Anno.csv") # manually annotated

plot_GREAT_anno(great2, n, "ARCHE2", anno)
plot_GREAT_anno(great5, n, "ARCHE5", anno)