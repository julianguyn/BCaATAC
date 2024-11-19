# Identify peaks significantly associated with drug response (univariable)

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(survcomp)
    library(ComplexHeatmap)
})


# load in cell line drug response data
load("DrugResponse/results/data/sensitivity_data.RData")

# load in and extract needed cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

# load in binary peak matrix
peak_mat <- fread("Signatures/data/bcacell_lines.tsv") |> as.data.frame()
rownames(peak_mat) <- paste(peak_mat$V1, peak_mat$V2, peak_mat$V3, sep = ":")
peak_mat <- peak_mat[,-c(1:3)] #num peaks: 1224557

# keep only needed samples and rename sample names
colnames(peak_mat) <- gsub("-", "\\.", colnames(peak_mat))
peak_mat <- peak_mat[,which(colnames(peak_mat) %in% samples$file)]
colnames(peak_mat) <- samples[match(colnames(peak_mat), samples$file),]$sample


### ===== Compute concordance index to find univariable associations ===== ###

# function to compute concordance index
computeCI <- function(signature_scores, sensitivity_data, label) {

    # keep and order common cell lines
    common_cells <- intersect(colnames(signature_scores), colnames(sensitivity_data))
    signature_scores <- signature_scores[,common_cells]
    sensitivity_data <- sensitivity_data[,common_cells]

    # create data frame to hold results
    combinations <- expand.grid(peak = rownames(signature_scores), drug = rownames(sensitivity_data))
    num_combinations <- nrow(combinations)

    # initialize columns to store univariable results
    combinations$ci <- NA
    combinations$pvalue <- NA
    combinations$se <- NA
    combinations$upper <- NA
    combinations$lower <- NA

    # extract drug and feature data into lists
    sen_list <- lapply(rownames(sensitivity_data), function(drug) as.numeric(sensitivity_data[drug, ]))
    feat_list <- lapply(rownames(signature_scores), function(peak) as.numeric(signature_scores[peak, ]))

    # compute concordance index for each feature-drug pair
    for (i in 1:num_combinations) {
        if (i %% 10000 == 0) { print(paste0(i, " out of ", num_combinations, " complete"))}
        
        # compute concordance index
        ci <- survcomp::concordance.index(
            sen_list[[which(rownames(sensitivity_data) == combinations$drug[i])]],                  # sensitivity data for all samples
            surv.time = feat_list[[which(rownames(signature_scores) == combinations$peak[i])]],     # feature values for all samples
            surv.event = rep(1, num_samples),
            outx = TRUE,
            method = "noether",
            na.rm = TRUE
        )

        combinations$ci[i] <- ci$c.index
        combinations$pvalue[i] <- ci$p.value
        combinations$se[i] <- ci$se
        combinations$upper[i] <- ci$upper
        combinations$lower[i] <- ci$lower
    }

    # filtering and multiple test correction
    combinations <- combinations[complete.cases(combinations$pvalue),]
    combinations$FDR <- p.adjust(combinations$pvalue, method = "BH")
    combinations$FDRsig <- combinations$FDR < 0.05

    # format dataframe for plotting (order by ci and add variable rank)
    combinations <- combinations[order(combinations$ci),]
    combinations$rank <- seq_len(nrow(combinations))
    combinations$pairs <- paste0(combinations$peak, "-", combinations$drug)
    combinations$pset <- label

    return(combinations)
}

# subset to compute only paclitaxel
ubr1_sen <- ubr1_sen["Paclitaxel",]
ubr2_sen <- ubr2_sen["Paclitaxel",]
gray_sen <- gray_sen["Paclitaxel",]
gcsi_sen <- gcsi_sen["Paclitaxel",]
gdsc_sen <- gdsc_sen["Paclitazel",]
ctrp_sen <- ctrp_sen["Paclitaxel",]
ccle_sen <- ccle_sen["Paclitaxel",]

# compute associations
ubr1_com <- computeCI(peak_mat, ubr1_sen, "UHNBreast1")
ubr2_com <- computeCI(peak_mat, ubr2_sen, "UHNBreast2")
gray_com <- computeCI(peak_mat, gray_sen, "GRAY")
gcsi_com <- computeCI(peak_mat, gcsi_sen, "gCSI")
gdsc_com <- computeCI(peak_mat, gdsc_sen, "GDSC")
ctrp_com <- computeCI(peak_mat, ctrp_sen, "CTRP")
ccle_com <- computeCI(peak_mat, ccle_sen, "CCLE")

# save results
save(ubr1_com, file = "../results/res1.RData")
save(ubr2_com, file = "../results/res2.RData")
save(gray_com, file = "../results/res3.RData")
save(gcsi_com, file = "../results/res4.RData")
save(gdsc_com, file = "../results/res5.RData")
save(ctrp_com, file = "../results/res6.RData")
save(ccle_com, file = "../results/res7.RData")

# save(ubr1_com, ubr2_com, gray_com, gcsi_com, gdsc_com, ctrp_com, ccle_com, file = "SupervisedDrugResponseAssociations.RData")


### ===== Assess concordance across PSets ===== ###

# function to filter for peaks with CI and FDR threshold
# CI threshold: CI > 0.8 or CI < 0.2
# FDR threshold: FDR < 0.05
filter_peaks <- function(combinations) {
    combinations <- combinations[which((combinations$ci > 0.7 | combinations$ci < 0.3) & combinations$FDR < 0.05),]
    return(combinations)
}

ubr1_com <- filter_peaks(ubr1_com) #86512
ubr2_com <- filter_peaks(ubr2_com) #137177
gray_com <- filter_peaks(gray_com) #149838
gcsi_com <- filter_peaks(gcsi_com) #68110
gdsc_com <- filter_peaks(gdsc_com) #61316
ctrp_com <- filter_peaks(ctrp_com) #90075
ccle_com <- filter_peaks(ccle_com) #71678

# upset plot of common peaks
set.seed(123)
toPlot <- make_comb_mat(list(
    UBR1 = ubr1_com$peak,
    UBR2 = ubr2_com$peak,
    GRAY = gray_com$peak,
    gCSI = gcsi_com$peak,
    GDSC = gdsc_com$peak,
    CTRP = ctrp_com$peak,
    CCLE = ccle_com$peak
))

# remove combinations of less than 10 pairs
toPlot <- toPlot[comb_size(toPlot) >= 1000]

# upset plot
pdf("SupervisedSignatures/results/figures/upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = c("UBR1", "UBR2", "GRAY", "gCSI", "GDSC", "CTRP", "CCLE"),
    comb_order = order(comb_size(toPlot), decreasing = TRUE),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()


### ===== Keep peaks with assocciations across PSets ===== ###

# merge dataframes
CI_combined <- rbind(ubr1_com, ubr2_com, gray_com, gcsi_com, gdsc_com, ctrp_com, ccle_com)

# keep peaks present in at least 4 PSets
peak_counts <- table(CI_combined$peak)
keep <- names(peak_counts[peak_counts > 4])
CI_combined <- CI_combined[CI_combined$peak %in% keep,]

# upset plot of common peaks
set.seed(123)
toPlot <- make_comb_mat(list(
    UBR1 = CI_combined$peak[CI_combined$pset == "UHNBreast1"],  #1774
    UBR2 = CI_combined$peak[CI_combined$pset == "UHNBreast2"],  #3003
    GRAY = CI_combined$peak[CI_combined$pset == "GRAY"],        #2840
    gCSI = CI_combined$peak[CI_combined$pset == "gCSI"],        #2655
    GDSC = CI_combined$peak[CI_combined$pset == "GDSC"],        #1127
    CTRP = CI_combined$peak[CI_combined$pset == "CTRP"],        #2169
    CCLE = CI_combined$peak[CI_combined$pset == "CCLE"]         #2336
))

# upset plot
pdf("SupervisedSignatures/results/figures/upset_feat_red.pdf", width=10, height=5)
UpSet(toPlot, set_order = c("UBR1", "UBR2", "GRAY", "gCSI", "GDSC", "CTRP", "CCLE"),
    comb_order = order(comb_size(toPlot), decreasing = TRUE),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()


# filter peak matrix to only keep peaks 
peak_mat <- peak_mat[rownames(peak_mat) %in% keep, ]
#save(peak_mat, file = "SupervisedSignatures/results/data/peak-filtered-peakmat.RData")


### ===== Remove correlated peaks ===== ###

# compute correlation matrix
peak_mat <- t(peak_mat)
corr_mat <- cor(peak_mat, use = "pairwise.complete.obs")
#save(corr_mat, file = "SupervisedSignatures/results/data/corr-mat.RData")

# set threshold
thres <- 0.8

# identify correlated peaks
upper_tri <- upper.tri(corr_mat)
correlated <- which(upper_tri & abs(corr_mat) > thres, arr.ind = TRUE)

# remove correlated peak from identified pairs
peak_reduced <- peak_mat[,-unique(correlated[,2])]


### ===== Prepare data for model training and testing ===== ###

# load in drug response data
load("SupervisedSignatures/results/data/pac-binarized.RData")

# remove cells with <NA> overall response
pac <- pac[!is.na(pac$response), ]
peak_reduced <- peak_reduced[rownames(peak_reduced) %in% rownames(pac),]

# order cell lines
pac <- pac[order(rownames(pac)),]
peak_reduced <- peak_reduced[order(rownames(peak_reduced)),]

write.csv(peak_reduced, file = "SupervisedSignatures/results/data/bca_peakmat.csv", quote = FALSE, row.names = TRUE)
write.csv(pac, file = "SupervisedSignatures/results/data/pac-binarized.csv", quote = FALSE, row.names = TRUE)