# Identify peaks significantly associated with drug response (univariable)

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(survcomp)
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
peak_mat <- peak_mat[,-c(1:3)]

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
    feat_list <- lapply(rownames(signature_scores), function(peak) as.numeric(-signature_scores[peak, ]))
    num_samples <- length(sen_list[[1]])

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
ubr2_sen <- ubr2_sen["Paclitaxel",]
gray_sen <- gray_sen["Paclitaxel",]
ctrp_sen <- ctrp_sen["Paclitaxel",]

# function to binarize drug repsonse
binarize_dr <- function(sensitivity_data) {
    return(as.data.frame(ifelse(sensitivity_data >= 0.5, 1, 0)))
}

# binarize drug response data
ubr1_sen <- binarize_dr(ubr1_sen)
ubr2_sen <- binarize_dr(ubr2_sen)
gray_sen <- binarize_dr(gray_sen)
gcsi_sen <- binarize_dr(gcsi_sen)
gdsc_sen <- binarize_dr(gdsc_sen)
ctrp_sen <- binarize_dr(ctrp_sen)
ccle_sen <- binarize_dr(ccle_sen)

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



### ===== Remove correlated peaks ===== ###

# on H4H: /cluster/projects/bhklab/projects/BCaATAC/SupervisedSignatures/procdata/peak-counts

# compute correlation matrix
peak_mat <- t(peak_mat)
dim(peak_mat) # 49     1224557
corr_mat <- cor(peak_mat, use = "pairwise.complete.obs")        # current dim: # 49     1224557
save(corr_mat, file = "corr-mat.RData")

upper_tri <- upper.tri(corr_mat)

# set threshold
thres <- 0.8
correlated <- which(upper_tri & abs(corr_mat) > thres, arr.ind = TRUE)

# remove correlated peak from identified pairs
peak_reduced <- peak_mat[-unique(correlated[,2]), ]

# Check the new dimensions of the matrix
dim(peak_reduced)

save(corr_mat, file = "corr-mat.RData")

