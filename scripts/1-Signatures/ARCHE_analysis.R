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
    library(readxl)
    library(ggsignif)
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

# load in metadata
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
meta$Sample.Name <- gsub("\\.", "-", meta$Sample.Name)

# get signature scores from NMF
mat <- get_arche_tcga()
mat$variable <- meta$Sample.Name[match(mat$variable, meta$ATAC.Seq.File.Name)]

# read in tcga clinical data
pheno <- read.table("data/rawdata/tcga/Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", header = TRUE, row.names = 1)
pheno <- t(pheno) |> as.data.frame()
rownames(pheno) <- gsub("\\.", "-", rownames(pheno))

# load in samples with rna and methylation
have_methylation <- readRDS("data/procdata/TCGA/samples_w_methylation.RDS")
have_rna <- readRDS("data/procdata/TCGA/samples_w_rna.RDS")
have_snv <- meta$Sample.Name[!is.na(meta$SNV.File.Name)]
dup_atac <- "TCGA-A2-A0T4"

###########################################################
# Plot ATAC-Signature heatmap
###########################################################

# add variables to matrix for plotting
mat$met <- ifelse(mat$variable %in% have_methylation, "Yes", "No")
mat$rna <- ifelse(mat$variable %in% have_rna, "Yes", "No")
mat$snv <- ifelse(mat$variable %in% have_snv, "Yes", "No")
mat$dup <- ifelse(mat$variable %in% dup_atac, "Yes", "No")

mat$age <- pheno$years_to_birth[match(mat$variable, rownames(pheno))] |> as.numeric()
mat$t_purity <- pheno$Tumor_purity[match(mat$variable, rownames(pheno))] |> as.numeric()
mat$stage <- pheno$pathologic_stage[match(mat$variable, rownames(pheno))]
mat$stageT <- pheno$pathology_T_stage[match(mat$variable, rownames(pheno))]
mat$stageN <- pheno$pathology_N_stage[match(mat$variable, rownames(pheno))]
mat$stageM <- pheno$pathology_M_stage[match(mat$variable, rownames(pheno))]
mat$OS <- pheno$overall_survival[match(mat$variable, rownames(pheno))] |> as.numeric()

# plot heatmap
mat$rank <- rep(1:75, each = 6)
plot_ARCHE_heatmap(mat)

###########################################################
# Plot staging variables
###########################################################

toPlot <- mat[,c("signature_assign", "stage", "stageT", "stageN", "stageM", "age", "OS", "t_purity")] |> unique()

# plot bar plots
plot_stage(toPlot, "stageT")
plot_stage(toPlot, "stageN")
plot_stage(toPlot, "stageM")
plot_stage(toPlot, "stage")

###########################################################
# ANOVA of clinical variables ~ ARCHE, and Tukey's HDS
###########################################################

aov(age ~ signature_assign, data = toPlot) |> summary()
#                 Df Sum Sq Mean Sq F value Pr(>F)
#signature_assign  5   1163   232.6   1.759  0.133
#Residuals        67   8863   132.3
aov(OS ~ signature_assign, data = toPlot) |> summary()
#                 Df   Sum Sq Mean Sq F value Pr(>F)
#signature_assign  5  2507859  501572   0.602  0.699
#Residuals        64 53324737  833199
tp.aov <- aov(t_purity ~ signature_assign, data = toPlot)
summary(tp.aov)
#                 Df Sum Sq Mean Sq F value  Pr(>F)
#signature_assign  5 0.3009 0.06018   6.447 5.7e-05 ***
#Residuals        69 0.6441 0.00933
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Tukey's HSD
tukey <- TukeyHSD(tp.aov)$signature_assign |> as.data.frame()
sig_tukey <- tukey[tukey$"p adj" < 0.05,]

###########################################################
# Plot other clinical variables
###########################################################

sig_codes <- round(sig_tukey$"p adj", 4)

plot_clinical(toPlot, "age")
plot_clinical(toPlot, "OS")
plot_clinical(toPlot, "t_purity")

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

###########################################################
# Plot HOMER findMotifsGenome results
###########################################################

plot_motifs("ARCHE1", "known")
plot_motifs("ARCHE1", "denovo")
plot_motifs("ARCHE2", "known")
plot_motifs("ARCHE2", "denovo")
plot_motifs("ARCHE3", "known")
plot_motifs("ARCHE3", "denovo")
plot_motifs("ARCHE4", "known")
plot_motifs("ARCHE4", "denovo")
plot_motifs("ARCHE5", "known")
plot_motifs("ARCHE5", "denovo")
plot_motifs("ARCHE6", "known")
plot_motifs("ARCHE6", "denovo")