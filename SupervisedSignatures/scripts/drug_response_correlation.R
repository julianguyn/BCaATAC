# Script to identify drugs of interest for supervised signature extraction

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
    library(plyr)
    library(ComplexHeatmap)
    library(circlize)
})

### ===== Get drugs of interest ===== ###

# read in PDX drug response data
pdx_response <- as.data.frame(read_excel("DrugResponsePDX/data/drugresponse/DrugResponse_PDX.xlsx", sheet = 1))
pdx_response[pdx_response == "NA"] <- NA
pdx_response[pdx_response == "TAXOL"] <- "PACLITAXEL"
pdx_response[pdx_response == "CARBOPLATIN"] <- "CARBOPLATINUM"
pdx_response <- na.omit(pdx_response)

# remove CONTROLs and H2O samples (left with n=23 drugs)
drugs_of_interest <- unique(pdx_response$drug[-which(pdx_response$drug %in% c("H2O", pdx_response$drug[grep("CONTROL", pdx_response$drug)]))])
pdx_response <- pdx_response[pdx_response$drug %in% drugs_of_interest,]

# keep only drugs with at least 10 samples (left with n=11 drugs)
drugs_of_interest <- names(table(pdx_response$drug)[table(pdx_response$drug) > 10])
pdx_response <- pdx_response[pdx_response$drug %in% drugs_of_interest,]


### ===== Get cell line drug response for drugs of interest ===== ### 

# load in cell line drug response data
load("DrugResponse/results/data/sensitivity_data.RData")

# load in and extract needed cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

# keep only drugs of interest
ubr1_sen <- ubr1_sen[toupper(rownames(ubr1_sen)) %in% drugs_of_interest,colnames(ubr1_sen) %in% samples$sample] #taxol only, 43 cells
ubr2_sen <- ubr2_sen[toupper(rownames(ubr2_sen)) %in% drugs_of_interest,colnames(ubr2_sen) %in% samples$sample] #taxol, carboplatin, & eribulin, 42 cells
gray_sen <- gray_sen[toupper(rownames(gray_sen)) %in% drugs_of_interest,colnames(gray_sen) %in% samples$sample] #taxol & carboplatin, 41 cells
gcsi_sen <- gcsi_sen[toupper(rownames(gcsi_sen)) %in% drugs_of_interest,colnames(gcsi_sen) %in% samples$sample] #taxol only, 25 cells
gdsc_sen <- gdsc_sen[toupper(rownames(gdsc_sen)) %in% drugs_of_interest,colnames(gdsc_sen) %in% samples$sample] #taxol only, 34 cells
ctrp_sen <- ctrp_sen[toupper(rownames(ctrp_sen)) %in% drugs_of_interest,colnames(ctrp_sen) %in% samples$sample] #taxol & carboplatin, 32 cells
ccle_sen <- ccle_sen[toupper(rownames(ccle_sen)) %in% drugs_of_interest,colnames(ccle_sen) %in% samples$sample] #taxol only, 23 cells


### ===== Correlate drug response for common cell lines ===== ### 

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

# get all paclitaxel response
pac <- rbind.fill(
    ubr1_sen[rownames(ubr1_sen) == "Paclitaxel",],
    ubr2_sen[rownames(ubr2_sen) == "Paclitaxel",],
    gray_sen[rownames(gray_sen) == "Paclitaxel",],
    gcsi_sen[rownames(gcsi_sen) == "Paclitaxel",],
    gdsc_sen[rownames(gdsc_sen) == "Paclitaxel",],
    ccle_sen[rownames(ccle_sen) == "Paclitaxel",],
    ctrp_sen[rownames(ctrp_sen) == "Paclitaxel",]
) |> t()
colnames(pac) <- c("UBR1", "UBR2", "GRAY", "gCSI", "GDSC", "CCLE", "CTRP")

# create correlation matrix
corr_mat <- cor(pac, use = "pairwise.complete.obs", method = "pearson")

# correlation heatmap for paclitaxel
png("SupervisedSignatures/results/figures/dr_corr/pcc-paclitaxel.png", width=125, height=100, units='mm', res = 600, pointsize=80)
plot_heatmap(corr_mat, "Pearson's")
dev.off()


# get all carboplatin response
car <- rbind.fill(
    ubr2_sen[rownames(ubr2_sen) == "Carboplatinum",],
    gray_sen[rownames(gray_sen) == "Carboplatinum",],
    ctrp_sen[rownames(ctrp_sen) == "Carboplatinum",]
) |> t()
colnames(car) <- c("UBR2", "GRAY", "CTRP")

# create correlation matrix
corr_mat <- cor(car, use = "pairwise.complete.obs", method = "pearson")

# correlation heatmap for carboplatinum
png("SupervisedSignatures/results/figures/dr_corr/pcc-carboplatinum.png", width=85, height=65, units='mm', res = 600, pointsize=80)
plot_heatmap(corr_mat, "Pearson's")
dev.off()


### ===== Concordance of binarized drug response values in cell lines ===== ### 

# function to binarize drug repsonse
binarize_dr <- function(sensitivity_data) {
    return(as.data.frame(ifelse(sensitivity_data >= 0.5, 1, 0)))
}

# binarize drug repsonse matrices
pac <- binarize_dr(pac)
car <- binarize_dr(car)


# function to compute jaccard similarity from: https://www.r-bloggers.com/2021/11/how-to-calculate-jaccard-similarity-in-r/
jaccard <- function(a, b) {
    # remove cell lines with NAs
    keep <- !is.na(a) & !is.na(b)
    a <- a[keep]
    b <- b[keep]

    intersection = sum(a == b)
    union = length(a)
    return (intersection/union)
}

# function to compute a matrix of jaccard similarity
jaccard_mat <- function(drug_response_mat) {

    # initialize matrix
    mat <- matrix(NA, nrow = ncol(drug_response_mat), ncol = ncol(drug_response_mat), 
                dimnames = list(names(drug_response_mat), names(drug_response_mat)))

    # compute pairwise jaccard similarity
    combn(seq_len(ncol(drug_response_mat)), 2, 
        function(idx) {
            i <- idx[1]
            j <- idx[2]
            
            # compute jaccard and fill in symmetric matrix
            mat[i, j] <<- mat[j, i] <<- jaccard(drug_response_mat[[i]], drug_response_mat[[j]])
        }
    )

    # set diagonal as 1
    diag(mat) <- 1

    return(mat)
}


# create jaccarad matrices
jacc_pac <- jaccard_mat(pac)
jacc_car <- jaccard_mat(car)

# correlation heatmap for paclitaxel
png("SupervisedSignatures/results/figures/dr_corr/jacc-paclitaxel.png", width=125, height=100, units='mm', res = 600, pointsize=80)
plot_heatmap(jacc_pac,"Jaccard\nSimilarity")
dev.off()

# correlation heatmap for paclitaxel
png("SupervisedSignatures/results/figures/dr_corr/jacc-carboplatinum.png", width=85, height=65, units='mm', res = 600, pointsize=80)
plot_heatmap(jacc_car,"Jaccard\nSimilarity")
dev.off()


# save binarized paclitaxel drug response
save(pac, file = "SupervisedSignatures/results/data/pac-binarized.RData")