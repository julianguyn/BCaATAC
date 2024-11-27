# Script to analyse the results from model training and testing on individual psets

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ComplexHeatmap)
})

set.seed(123)

# function to read in individual feature importance dataframes
read_features <- function(path) {

    # load in files
    files <- list.files(path, recursive = T, pattern = "features")

    # create data frame to hold results
    features <- as.data.frame(matrix(nrow = 0, ncol = 4))
    colnames(features) <- c("Peak", "Weight", "Model", "PSet")

    for (file in files) {

        # extract model and pset
        group <- gsub("/.*", "", file)  # will either be model (indiv pset) or 'averaged'

        if (group == "averaged") {
            pset = group
            model = gsub("averaged/", "", gsub("_.*", "", file))
        } else {
            model = group
            pset <- unlist(strsplit(file, "_"))[2]
        }

        # read in feature importance file
        df <- read.csv(paste0(path, file))

        # keep only peaks with a non-zero weight
        df <- df[df$Weight != 0,]

        # include model and pset details
        df$Model <- model
        df$PSet <- pset

        # save results
        features <- rbind(features, df)
    }

    return(features)
}

# function to threshold peaks by model weights and plot distribution
thres_weight <- function(features, model, thres) {

    # subset by model of interest
    df <- features[features$Model == model,]

    # format df for plotting
    df <- df[order(df$Weight, decreasing = T),]
    df$Rank <- 1:nrow(df)

    # annotate points with weight greater than threshold
    df$Sig <- ifelse(abs(df$Weight) >= thres, 'Above Thres', 'Below Thres')

    # get distribution of number of peaks
    table(df$PSet, df$Sig) |> print()

    # plot distribution of weights
    p <- ggplot(df) + geom_point(aes(x = Rank, y = Weight, color = Sig)) +
        facet_wrap(~PSet, nrow = 2) + geom_hline(yintercept = 0, linetype = 'dotted') +
        scale_color_manual(values = c("#7294D4", "#BC4749"), limits = c('Above Thres', 'Below Thres')) +
        theme_classic() + theme(
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            axis.ticks.x = element_blank(), 
            axis.text.x = element_blank(),
            legend.key.size = unit(0.7, 'cm')
        ) + labs(x = "Peak", color = "")

    return(p)
}

# function to keep only top n peaks from each PSet
top_features <- function(features, n) {

    # get unique groups
    features$label <- paste0(features$Model, "-", features$PSet)
    groups <- unique(features$label)

    # create data frame to hold results
    top_features <- as.data.frame(matrix(nrow = 0, ncol = 4))
    colnames(top_features) <- c("Peak", "Weight", "Model", "PSet")
    
    # loop through each pset
    for (group in groups) {

        # subset for pset and order features
        sub <- features[features$label == group,]
        sub <- sub[order(abs(sub$Weight), decreasing = T),]

        # if >50 features, keep top 50; otherwise keep all
        if (nrow(sub) > n) { 
            sub <- sub[1:n,] 
        } 

        # save features
        top_features <- rbind(top_features, sub)
    }

    return(top_features)
}


### ===== AAC Drug Response - Top 50 Peaks ===== ###

# read in feature importance from AAC
AAC_path <- "SupervisedSignatures/results/data/AAC/"
AAC_features <- read_features(AAC_path)

# keep top 50 peaks from each PSet
AAC_features <- top_features(AAC_features, 50)
AAC_features$label <- paste0(AAC_features$Model, "-", AAC_features$PSet)

# create upset plots of overlapping 'sig' peaks
toPlot <- make_comb_mat(list(
    LM_UBR1 = AAC_features$Peak[AAC_features$label == "lm-UBR1"],   LM_UBR2 = AAC_features$Peak[AAC_features$label == "lm-UBR2"],
    LM_GRAY = AAC_features$Peak[AAC_features$label == "lm-GRAY"],   LM_GDSC = AAC_features$Peak[AAC_features$label == "lm-GDSC"],
    LM_CTRP = AAC_features$Peak[AAC_features$label == "lm-CTRP"],   LM_CCLE = AAC_features$Peak[AAC_features$label == "lm-CCLE"],
    LM_AVG = AAC_features$Peak[AAC_features$label == "lm-averaged"] 
))
pdf("SupervisedSignatures/results/figures/indiv_top/AAC_lm_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()

toPlot <- make_comb_mat(list(
    Lasso_GRAY = AAC_features$Peak[AAC_features$label == "lasso-GRAY"],     Lasso_gCSI = AAC_features$Peak[AAC_features$label == "lasso-gCSI"],
    Lasso_GDSC = AAC_features$Peak[AAC_features$label == "lasso-GDSC"],     Lasso_CTRP = AAC_features$Peak[AAC_features$label == "lasso-CTRP"],   
    Lasso_CCLE = AAC_features$Peak[AAC_features$label == "lasso-CCLE"],     Lasso_AVG = AAC_features$Peak[AAC_features$label == "lasso-averaged"]
)) 
pdf("SupervisedSignatures/results/figures/indiv_top/AAC_lasso_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()

toPlot <- make_comb_mat(list(  
    Ridge_UBR2 = AAC_features$Peak[AAC_features$label == "ridge-UBR2"], Ridge_GRAY = AAC_features$Peak[AAC_features$label == "ridge-GRAY"],
    Ridge_gCSI = AAC_features$Peak[AAC_features$label == "ridge-gCSI"], Ridge_CTRP = AAC_features$Peak[AAC_features$label == "ridge-CTRP"],
    Ridge_CCLE = AAC_features$Peak[AAC_features$label == "ridge-CCLE"], Ridge_AVG = AAC_features$Peak[AAC_features$label == "ridge-averaged"]
))
pdf("SupervisedSignatures/results/figures/indiv_top/AAC_ridge_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()

toPlot <- make_comb_mat(list(
    EN_UBR1 = AAC_features$Peak[AAC_features$label == "en-UBR1"],   EN_UBR2 = AAC_features$Peak[AAC_features$label == "en-UBR2"],
    EN_GRAY = AAC_features$Peak[AAC_features$label == "en-GRAY"],   EN_gCSI = AAC_features$Peak[AAC_features$label == "en-gCSI"],
    EN_CCLE = AAC_features$Peak[AAC_features$label == "en-CCLE"],   EN_AVG = AAC_features$Peak[AAC_features$label == "en-averaged"]
))
pdf("SupervisedSignatures/results/figures/indiv_top/AAC_en_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()

toPlot <- make_comb_mat(list(
    RF_UBR1 = AAC_features$Peak[AAC_features$label == "rf-UBR1"],   RF_UBR2 = AAC_features$Peak[AAC_features$label == "rf-UBR2"],
    RF_GRAY = AAC_features$Peak[AAC_features$label == "rf-GRAY"],   RF_gCSI = AAC_features$Peak[AAC_features$label == "rf-gCSI"],
    RF_GDSC = AAC_features$Peak[AAC_features$label == "rf-GDSC"],   RF_CTRP = AAC_features$Peak[AAC_features$label == "rf-CTRP"],  
    RF_CCLE = AAC_features$Peak[AAC_features$label == "rf-CCLE"],   RF_AVG = AAC_features$Peak[AAC_features$label == "rf-averaged"]
))
pdf("SupervisedSignatures/results/figures/indiv_top/AAC_rf_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()


### ===== Binarized Drug Response - Top 50 Peaks ===== ###

# read in feature importance from binarized drug response
bin_path <- "SupervisedSignatures/results/data/bin/"
bin_features <- read_features(bin_path)

# keep top 50 peaks from each PSet
bin_features <- top_features(bin_features, 50)
bin_features$label <- paste0(bin_features$Model, "-", bin_features$PSet)

# create upset plot of overlapping 'sig' peaks
toPlot <- make_comb_mat(list(
    Log_UBR1 = bin_features$Peak[bin_features$label == "log-UBR1"], Log_UBR2 = bin_features$Peak[bin_features$label == "log-UBR2"],
    Log_GRAY = bin_features$Peak[bin_features$label == "log-GRAY"], Log_gCSI = bin_features$Peak[bin_features$label == "log-gCSI"],
    Log_CTRP = bin_features$Peak[bin_features$label == "log-CTRP"], Log_AVG = bin_features$Peak[bin_features$label == "log-averaged"]
))
pdf("SupervisedSignatures/results/figures/indiv_top/bin_log_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()

toPlot <- make_comb_mat(list(
    Lasso_UBR1 = bin_features$Peak[bin_features$label == "lasso-UBR1"], Lasso_GRAY = bin_features$Peak[bin_features$label == "lasso-GRAY"],
    Lasso_gCSI = bin_features$Peak[bin_features$label == "lasso-gCSI"], Lasso_CTRP = bin_features$Peak[bin_features$label == "lasso-CTRP"],
    Lasso_AVG = bin_features$Peak[bin_features$label == "lasso-averaged"]
))
pdf("SupervisedSignatures/results/figures/indiv_top/bin_lasso_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()

toPlot <- make_comb_mat(list(   
    Ridge_UBR1 = bin_features$Peak[bin_features$label == "ridge-UBR1"], Ridge_UBR2 = bin_features$Peak[bin_features$label == "ridge-UBR2"],
    Ridge_GRAY = bin_features$Peak[bin_features$label == "ridge-GRAY"], Ridge_gCSI = bin_features$Peak[bin_features$label == "ridge-gCSI"],
    Ridge_CTRP = bin_features$Peak[bin_features$label == "ridge-CTRP"], Ridge_AVG = bin_features$Peak[bin_features$label == "ridge-averaged"]
))
pdf("SupervisedSignatures/results/figures/indiv_top/bin_ridge_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()

toPlot <- make_comb_mat(list(
    EN_UBR1 = bin_features$Peak[bin_features$label == "en-UBR1"],   EN_UBR2 = bin_features$Peak[bin_features$label == "en-UBR2"],
    EN_GRAY = bin_features$Peak[bin_features$label == "en-GRAY"],   EN_CTRP = bin_features$Peak[bin_features$label == "en-CTRP"],
    EN_AVG = bin_features$Peak[bin_features$label == "en-averaged"]
))
pdf("SupervisedSignatures/results/figures/indiv_top/bin_en_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()

toPlot <- make_comb_mat(list(
    RF_UBR1 = bin_features$Peak[bin_features$label == "rf-UBR1"],   RF_UBR2 = bin_features$Peak[bin_features$label == "rf-UBR2"],
    RF_GRAY = bin_features$Peak[bin_features$label == "rf-GRAY"],   RF_gCSI = bin_features$Peak[bin_features$label == "rf-gCSI"],
    RF_CTRP = bin_features$Peak[bin_features$label == "rf-CTRP"],   RF_AVG = bin_features$Peak[bin_features$label == "rf-averaged"]
))
pdf("SupervisedSignatures/results/figures/indiv_top/bin_rf_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()



# TODO: update the following section (if we want to keep it?)

### ===== AAC Drug Response - Thresholding ===== ###

# read in feature importance from AAC
AAC_path <- "SupervisedSignatures/results/data/AAC/"
AAC_features <- read_features(AAC_path)

# find and plot thresholds for peak weights
png("SupervisedSignatures/results/figures/indiv_thres/AAC-lm.png", width=200, height=150, units='mm', res = 600, pointsize=80)
thres_weight(AAC_features, 'lm', thres = 0.005)
dev.off()

png("SupervisedSignatures/results/figures/indiv_thres/AAC-ridge.png", width=200, height=150, units='mm', res = 600, pointsize=80)
thres_weight(AAC_features, 'ridge', thres = 0.0025)
dev.off()

# subset feature dataframe by threshold
AAC_features <- AAC_features[-which(AAC_features$Model == 'lm' & abs(AAC_features$Weight) < 0.005),] 
AAC_features <- AAC_features[-which(AAC_features$Model == 'ridge' & abs(AAC_features$Weight) < 0.0025),]
AAC_features <- na.omit(AAC_features)


# create upset plot of overlapping 'sig' peaks
AAC_features$label <- paste0(AAC_features$Model, "-", AAC_features$PSet)
toPlot <- make_comb_mat(list(
    LM_UBR1 = AAC_features$Peak[AAC_features$label == "lm-UBR1"],   LM_UBR2 = AAC_features$Peak[AAC_features$label == "lm-UBR2"],
    LM_GRAY = AAC_features$Peak[AAC_features$label == "lm-GRAY"],   LM_CTRP = AAC_features$Peak[AAC_features$label == "lm-CTRP"],
    LM_CCLE = AAC_features$Peak[AAC_features$label == "lm-CCLE"],
    
    Ridge_UBR2 = AAC_features$Peak[AAC_features$label == "ridge-UBR2"], Ridge_GRAY = AAC_features$Peak[AAC_features$label == "ridge-GRAY"],
    Ridge_gCSI = AAC_features$Peak[AAC_features$label == "ridge-gCSI"], Ridge_CTRP = AAC_features$Peak[AAC_features$label == "ridge-CTRP"],
    Ridge_CCLE = AAC_features$Peak[AAC_features$label == "ridge-CCLE"],

    EN_UBR1 = AAC_features$Peak[AAC_features$label == "en-UBR1"],   EN_UBR2 = AAC_features$Peak[AAC_features$label == "en-UBR2"],
    EN_GRAY = AAC_features$Peak[AAC_features$label == "en-GRAY"],   EN_gCSI = AAC_features$Peak[AAC_features$label == "en-gCSI"],
    #EN_CTRP = AAC_features$Peak[AAC_features$label == "en-CTRP"],  
    EN_CCLE = AAC_features$Peak[AAC_features$label == "en-CCLE"]
))


# upset plot
pdf("SupervisedSignatures/results/figures/indiv_thres/AAC_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()



### ===== Binarized Drug Response - Thresholding ===== ###

# read in feature importance from AAC
bin_path <- "SupervisedSignatures/results/data/bin/"
bin_features <- read_features(bin_path)

# find and plot thresholds for peak weights
png("SupervisedSignatures/results/figures/indiv_thres/bin-log.png", width=200, height=150, units='mm', res = 600, pointsize=80)
thres_weight(bin_features, 'log', thres = 0.2)
dev.off()

png("SupervisedSignatures/results/figures/indiv_thres/bin-ridge.png", width=200, height=150, units='mm', res = 600, pointsize=80)
thres_weight(bin_features, 'ridge', thres = 0.05)
dev.off()

png("SupervisedSignatures/results/figures/indiv_thres/bin-en.png", width=200, height=150, units='mm', res = 600, pointsize=80)
thres_weight(bin_features, 'en', thres = 0.05)
dev.off()

# subset feature dataframe by threshold
bin_features <- bin_features[-which(bin_features$Model == 'log' & abs(bin_features$Weight) < 0.2),] 
bin_features <- bin_features[-which(bin_features$Model == 'ridge' & abs(bin_features$Weight) < 0.05),]
bin_features <- bin_features[-which(bin_features$Model == 'en' & abs(bin_features$Weight) < 0.05),] 
bin_features <- na.omit(bin_features)


# create upset plot of overlapping 'sig' peaks
bin_features$label <- paste0(bin_features$Model, "-", bin_features$PSet)
toPlot <- make_comb_mat(list(
    Log_UBR1 = bin_features$Peak[bin_features$label == "log-UBR1"], Log_UBR2 = bin_features$Peak[bin_features$label == "log-UBR2"],
    Log_GRAY = bin_features$Peak[bin_features$label == "log-GRAY"], Log_gCSI = bin_features$Peak[bin_features$label == "log-gCSI"],
    Log_CTRP = bin_features$Peak[bin_features$label == "log-CTRP"],

    Lasso_UBR1 = bin_features$Peak[bin_features$label == "lasso-UBR1"], Lasso_GRAY = bin_features$Peak[bin_features$label == "lasso-GRAY"],
    Lasso_gCSI = bin_features$Peak[bin_features$label == "lasso-gCSI"], Lasso_CTRP = bin_features$Peak[bin_features$label == "lasso-CTRP"],
    
    Ridge_UBR1 = bin_features$Peak[bin_features$label == "ridge-UBR1"], Ridge_UBR2 = bin_features$Peak[bin_features$label == "ridge-UBR2"],
    Ridge_GRAY = bin_features$Peak[bin_features$label == "ridge-GRAY"],

    EN_UBR1 = bin_features$Peak[bin_features$label == "en-UBR1"],   EN_UBR2 = bin_features$Peak[bin_features$label == "en-UBR2"],
    EN_GRAY = bin_features$Peak[bin_features$label == "en-GRAY"]
))


# upset plot
pdf("SupervisedSignatures/results/figures/indiv_thres/bin_features_upset.pdf", width=10, height=5)
UpSet(toPlot, set_order = rownames(toPlot),
    comb_order = order(comb_size(toPlot), decreasing = F),
    top_annotation = upset_top_annotation(toPlot, add_numbers = TRUE))
dev.off()



