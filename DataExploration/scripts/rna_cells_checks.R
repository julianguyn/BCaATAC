# Compare the RNA-Seq expression to other PSets

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(ggplot2)
    library(reshape2)
    library(ggpubr)
    library(ComplexHeatmap)
    library(circlize)
})

set.seed(101)

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
samples$file <- gsub("\\.", "-", gsub("\\.$", "", samples$file))
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

save(cells, gray, gcsi, ccle, file = "DataExploration/data/rnaseq_melted.RData")

### ===== Compute correlations between RNA-Seq data of PSets ===== ###

# function to create a correlation plot and compute correlation coefficient between two psets
corr_pset <- function(pset1, pset2, lab1, lab2) {

    # intersected cell-gene pairs
    pset1$pairs <- paste0(pset1$Var1, pset1$Var2)
    pset2$pairs <- paste0(pset2$Var1, pset2$Var2)
    
    overlap <- intersect(pset1$pairs, pset2$pairs)
    
    # save intersected drug-cell pairs and order them
    pset1_intersect <- pset1[pset1$pairs %in% overlap,]
    pset1_intersect <- pset1[match(overlap, pset1$pairs),]
    pset2_intersect <- pset2[pset2$pairs %in% overlap,]
    pset2_intersect <- pset2[match(overlap, pset2$pairs),]

    # pearson correlation coefficient
    df <- data.frame(pair = overlap, pset1 = pset1_intersect$value, pset2 = pset2_intersect$value)
    corr <- cor(df$pset1, df$pset2,  method = "pearson", use = "complete.obs")
    print(paste("Correlation between", lab1, "and", lab2, ":", corr))

    return(corr)
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

# initialize correlation matrix
psets <- c("UHNBreast2", "GRAY", "gCSI", "CCLE")
mat <- matrix(NA, nrow = 4, ncol = 4, dimnames = list(psets, psets))

# fill in matrix values
mat[1,2] <- mat[2,1] <- p1
mat[1,3] <- mat[3,1] <- p2
mat[1,4] <- mat[4,1] <- p3
mat[2,3] <- mat[3,2] <- p4
mat[2,4] <- mat[4,2] <- p5
mat[3,4] <- mat[4,3] <- p6

# set diagonal as 1
diag(mat) <- 1

save(mat, file = "DataExploration/data/corr_mat.RData")

# set palette for plotting
col_fun <- colorRamp2(c(-1, 0, 1), c("#A85751", "white", "#66999B"))

# function to create heatmap
plot_heatmap <- function(mat, legend_label) {
    plot <- Heatmap(mat,
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        heatmap_legend_param = list(
            title = legend_label,
            color_bar = "continuous"
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(round(mat[i, j], 2), x, y, gp = gpar(fontsize = 8))
        }
    )
    return(plot)
}

png("DataExploration/results/figures/pearsons_corr.png", width=125, height=100, units='mm', res = 600, pointsize=80)
plot_heatmap(mat, "Pearson's")
dev.off()


### ===== Identify Highly Correlated Cell Lines from UBR2 ===== ###

# function to correlate cell lines from UBR2 to cell lines in other psets
corr_cells <- function(pset, pset_label) {

    # get unique cell lines from UHNBreast2 and pset
    ubr2_cells <- unique(cells$Var2)
    pset_cells <- unique(pset$Var2)

    # initalize df to store results
    num_combo <- length(ubr2_cells) * length(pset_cells)
    corr_res <- data.frame(matrix(nrow=num_combo, ncol=3))
    colnames(corr_res) <- c("UBR2_Cell", "PSet_Cell", "Pearsons")

    # set iteration
    i = 0

    # loop through for every cell line in UBR2
    for (cl in ubr2_cells) {

        # subset for unique cell line
        subset_ubr2 <- cells[cells$Var2 == cl,]

        # loop through for every cell line in pset
        for (pset_cl in pset_cells) {

            # increment iteration
            i = i + 1
            #if (i %% 20 == 0) {print(paste0(i, " out of ", num_combo, " complete"))}

            # subset for this unique cell line
            subset_pset <- pset[pset$Var2 == pset_cl,]

            # check genes are ordered same TODO: print message if not all TRUE
            #print(table(subset_ubr2$Var1 == subset_pset$Var1))

            # compute pearson's correlation coefficient
            corr <- cor(subset_ubr2$value, subset_pset$value,  method = "pearson", use = "complete.obs")

            # store correlation results in dataframe
            corr_res$UBR2_Cell[i] <- cl
            corr_res$PSet_Cell[i] <- pset_cl
            corr_res$Pearsons[i] <- corr
        }
    }

    corr_res$PSet <- pset_label
    return(corr_res)
}


# rename replicate in gray
gray$Var2 <- as.character(gray$Var2)
gray$Var2[1380149:1429439] <- rep("UACC812_r2", 49291)
corr_gray <- corr_cells(gray, "GRAY")
corr_gcsi <- corr_cells(gcsi, "gCSI")
corr_ccle <- corr_cells(ccle, "CCLE")

save(corr_gray, corr_gcsi, corr_ccle, file = "DataExploration/data/rnaseq-corr.RData")


# function to identify top correlations (num_top) for each UBR2 cell line
top_corr <- function(corr_df, num_top) {

    # get unique cell lines from UHNBreast2
    ubr2_cells <- unique(corr_df$UBR2_Cell)

    # initalize df to store results
    top_res <- data.frame(matrix(NA, nrow = 0, ncol=4))
    colnames(top_res) <- colnames(corr_df)

    # loop through each unique cell line
    for (cl in ubr2_cells) {

        # subset for unique cell line
        subset_ubr2 <- corr_df[corr_df$UBR2_Cell == cl,]

        # order and extract top num_top pearson's correlation
        subset_ubr2 <- subset_ubr2[order(subset_ubr2$Pearsons, decreasing = T),]
        subset_ubr2 <- subset_ubr2[1:num_top,]

        # save results
        top_res <- rbind(top_res, subset_ubr2)
    }

    return(top_res)
}

# get top correlations
gray_top <- top_corr(corr_gray, num_top = 1)
gcsi_top <- top_corr(corr_gcsi, num_top = 1)
ccle_top <- top_corr(corr_ccle, num_top = 1)

# merge and melt results for plotting
merge_df <- rbind(gray_top, gcsi_top, ccle_top)
merge_df <- melt(merge_df)

# plot top cell line match
png("DataExploration/results/figures/cl_correlations.png", width=150, height=150, units='mm', res = 600, pointsize=80)
ggplot(merge_df, aes(x = PSet, y = UBR2_Cell, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = PSet_Cell), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "#008DD5", limits = c(0, 0.5)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(0.5, 'cm')) +
  labs(fill = "Pearsons", y = "UHNBreast2 Cell Line")
dev.off()