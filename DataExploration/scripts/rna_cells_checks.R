# Compare the RNA-Seq expression to other PSets

setwd("C:/Users/julia/Documents/BCaATAC")

library(PharmacoGx)
library(umap)
library(ggplot2)
library(reshape2)

set.seed(101)

# set up palette for plotting
subtype_pal <- c("Basal" = "#394032", "Her2" = "#A6A57A","LumA" = "#8C5E58","LumB" = "#5A352A","Normal" = "#8F8073", "#eFeBF7")

# load in cell line RNA seq counts
uhnbreast2 <- readRDS("DrugResponse/data/PharmacoSet.RDS")
cells <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@assays@data$expr
colnames(cells) <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@colData$sampleid

# load in cell line annotation
cells_meta <- read.csv("MetaData/Lupien/BCa_samples.csv")
# from https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2911-z/tables/3, set 600MPE to LuminalA
cells_meta[which(cells_meta$subtype == "cell_line"),]$subtype <- "LumA"

# keep only cell lines being used
samples <- read.csv("DrugResponse/data/cl_samples.csv")
samples$file <- gsub("\\.$", "", samples$file)
samples$file <- gsub("\\.", "-", samples$file)
samples$match <- gsub("-", "", samples$sample)

# remove dups
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

cells_meta <- cells_meta[which(cells_meta$sample %in% samples$file),]
cells <- cells[,colnames(cells) %in% samples$match]

save(cells, file = "DataExploration/data/cells_rnaseq.RData")

###########################
### Load in other PSETS ###
###########################

# load in Gray PSet
gray <- readRDS("DrugResponse/data/PSet_GRAY2017.rds")
gray <- updateObject(gray)

# get rnaseq counts
counts <- gray@molecularProfiles$Kallisto_0.46.1.rnaseq.counts@assays@data$exprs
colnames(counts) <- gsub("-", "", gray@molecularProfiles$Kallisto_0.46.1.rnaseq.counts@colData$sampleid)
gray <- counts[,colnames(counts) %in% colnames(cells)]
gray <- gray[rownames(gray) %in% rownames(cells),]


# load in gCSI PSet
gcsi <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/gCSI.rds")
gcsi <- updateObject(gcsi)

# get rnaseq counts
counts <- gcsi@molecularProfiles$Kallisto_0.46.1.rnaseq.counts@assays@data$exprs
colnames(counts) <- gsub("-", "", gcsi@molecularProfiles$Kallisto_0.46.1.rnaseq.counts@colData$sampleid)
gcsi <- counts[,colnames(counts) %in% colnames(cells)]
gcsi <- gcsi[rownames(gcsi) %in% rownames(cells),]


# load in GDSC PSet #only 4 overlap?

# load in CCLE PSet
ccle <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/CCLE.rds")
ccle <- updateObject(ccle)

# get rnaseq counts
counts <- ccle@molecularProfiles$Kallisto_0.46.1.rnaseq.counts@assays@data$exprs
colnames(counts) <- gsub("-", "", ccle@molecularProfiles$Kallisto_0.46.1.rnaseq.counts@colData$sampleid)
ccle <- counts[,colnames(counts) %in% colnames(cells)]
ccle <- ccle[rownames(ccle) %in% rownames(cells),]

# keep overlapping genes
cells <- cells[rownames(cells) %in% rownames(ccle),]

# format rnaseq counts for correlation
cells <- melt(cells)
gray <- melt(gray)
gcsi <- melt(gcsi)
ccle <- melt(ccle)


# function to create a correlation plot and compute correlation coefficient between two psets
corr_pset <- function(pset1, pset2, lab1, lab2) {
    ################################################
    # Inputs:
    #   pset1: rna counts of first pset
    #   pset2: rna counts of second pset
    #   lab1: label of first pset
    #   lab2: label of second pset
    # Outputs:
    #   p: correlation graph 
    #   corr_df: will update this dataframe with 
    #            Pearson's correlation coefficient
    #            and number of overlapping features
    ################################################

    # intersected cell-gene pairs
    pset1$pairs <- paste0(pset1$Var1, pset1$Var2)
    pset2$pairs <- paste0(pset2$Var1, pset2$Var2)
    
    overlap <- intersect(pset1$pairs, pset2$pairs)
    
    # save intersected drug-cell pairs and order them
    pset1_intersect <- pset1[pset1$pairs %in% overlap,]
    pset1_intersect <- pset1[match(overlap, pset1$pairs),]
    pset2_intersect <- pset2[pset2$pairs %in% overlap,]
    pset2_intersect <- pset2[match(overlap, pset2$pairs),]

    # scatter plot of drug response difference for CTRP and GDSC 
    toPlot <- data.frame(pair = overlap, pset1 = pset1_intersect$value, pset2 = pset2_intersect$value)
    
    # pearson correlation coefficient
    corr <- cor(toPlot$pset1, toPlot$pset2,  method = "pearson", use = "complete.obs")

    p <- ggplot(toPlot, aes(x = pset1, y = pset2)) + 
            geom_smooth(method=lm, show.legend = FALSE, color = "#046C9A") + geom_point(shape = 21, size = 2.5, color = "black", fill = "#899DA4") + 
            geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") + 
            theme_classic() + xlim(c(0, 1)) + ylim(c(0, 1)) + theme(legend.key.size = unit(0.5, 'cm')) +
            labs(x = lab1, y = lab2) + geom_text(x = 0.01, y = 1, label = paste("Corr: ", round(corr, digits = 3)), color = "black")

    print(paste0("Correlation between ", lab1, " and ", lab2, ": ", corr))

    return(p)
}

# UHN Breast2
p1 <- corr_pset(cells, gray, "UHNBreast2", "GRAY") #0.265757812435099
p2 <- corr_pset(cells, gcsi, "UHNBreast2", "gCSI") #0.261650434507992
p3 <- corr_pset(cells, ccle, "UHNBreast2", "CCLE") #0.252381917686653

# GRAY
p4 <- corr_pset(gray, gcsi, "GRAY", "gCSI") #0.918926917257591
p5 <- corr_pset(gray, ccle, "GRAY", "CCLE") #0.915262101138858

# gCSI
p6 <- corr_pset(gcsi, ccle, "gCSI", "CCLE") #0.943145494132632


suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

png("DataExploration/results/figures/rna_cells_checks.png", width = 15, height = 13, res = 600, units = "in")
grid.arrange(p1, p2, p3, p4, p5, p6,
             ncol = 3, nrow = 3,
             layout_matrix = rbind(c(1, NA, NA),
                                   c(2,  4, NA),
                                   c(3,  5, 6)))
dev.off()


cp HDQP-1_2_S25_L002_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp HS-578-T_2_S1_L001_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp HS578._peaks.narrowPeak /cluster/home/julian/temp/peaks
cp JIMT-1_2_S3_L000_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp KPL1_2_S15_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp LY2_1_S2_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp MB157._peaks.narrowPeak /cluster/home/julian/temp/peaks
cp MCF10A_1_S2_L001_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp MCF7_1_S14_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp MDA_MB_175_VII_2_S5_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp MDA-MB-231._peaks.narrowPeak /cluster/home/julian/temp/peaks
cp MDA-MB-436._peaks.narrowPeak /cluster/home/julian/temp/peaks
cp MDA-MB-468._peaks.narrowPeak /cluster/home/julian/temp/peaks
cp MDA231_1_S2_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp MDA361_2_S16_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp MDA468_2_S11_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp MPE600_2_S10_L001_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp MX1._peaks.narrowPeak /cluster/home/julian/temp/peaks
cp SKBR-7_2_S17_L002_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp SKBR3_2_S8_L001_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp SKBR5_2_S2_L001_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp SUM149_2_S15_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp SUM149PT._peaks.narrowPeak /cluster/home/julian/temp/peaks
cp SUM159PT_2_S8_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp SUM159PT._peaks.narrowPeak /cluster/home/julian/temp/peaks
cp SUM52PE_1_S7_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp T47D_2_S14_L002_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp UACC812_2_S14_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
cp ZR-75-1_2_S1_L001_peaks.filtered.narrowPeak /cluster/home/julian/temp/peaks
 
