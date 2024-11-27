# Plot model results across folds

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(reshape2)
    library(ggpubr)
})

# function to plot correlations (for AAC)
plot_folds <- function(df, avg_df, type) {              
    # type: 'AAC' or 'bin'

    # format avg_df and merge two dataframes for plotting
    avg_df$PSet <- "Averaged"
    df <- rbind(df, avg_df)

    # specify response type dependent variables
    if (type == 'AAC') {
        keep <- c("PSet", "Fold", "Spearman", "Pearson")
        pal = c("#603A40", "#84596B")
        yint = c(0.7, 0.8, 0.9)
        label = "Correlation"
    } else {
        keep <- c("PSet", "Fold", "Precision", "Recall", "F1")
        pal = c("#603A40", "#84596B", "#9A7E8F")
        yint = c(0.67, 0.8, 1)
        label = "Score"
    }

    # keep only necessary columns and replace NA with 0
    df <- df[,keep]
    df[is.na(df)] <- 0

    # format dataframe for plotting
    df$Fold <- paste(df$PSet, df$Fold, sep = "-")
    toPlot <- reshape2::melt(df)
    toPlot$Fold <- gsub(".*-", "", toPlot$Fold)

    # pset options
    psets <- c("UBR1", "UBR2", "GRAY", "gCSI", "GDSC", "CCLE", "CTRP", "Averaged")
    toPlot$PSet <- factor(toPlot$PSet, levels = psets[psets %in% unique(toPlot$PSet)])

    p <- ggplot(toPlot, aes(x = Fold, y = value, fill = variable)) + 
            geom_bar(stat = "identity", position = 'dodge') +
            facet_wrap(~PSet, nrow = 1) +   scale_fill_manual(values = pal) +
            theme_classic() + 
            geom_hline(yintercept = 0, linetype = "solid") + 
            geom_hline(yintercept = yint, linetype = "dotted") +
            theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
            labs(y = label, fill = label)
    
    return(p)
}

### ===== AAC Model Performance Evaluation across Folds ===== ###

AAC_path <- "SupervisedSignatures/results/data/AAC/"

# load in indiv pset results
lm <- read.csv(paste0(AAC_path, "lm/lm.csv"))
lasso <- read.csv(paste0(AAC_path, "lasso/lasso.csv"))
ridge <- read.csv(paste0(AAC_path, "ridge/ridge.csv"))
en <- read.csv(paste0(AAC_path, "en/en.csv"))
rf <- read.csv(paste0(AAC_path, "rf/rf.csv"))  

# load in averaged results
avg_lm <- read.csv(paste0(AAC_path, "averaged/lm.csv"))
avg_lasso <- read.csv(paste0(AAC_path, "averaged/lasso.csv"))
avg_ridge <- read.csv(paste0(AAC_path, "averaged/ridge.csv"))
avg_en <- read.csv(paste0(AAC_path, "averaged/en.csv"))
avg_rf <- read.csv(paste0(AAC_path, "averaged/rf.csv"))

# plot results across folds
p1 <- plot_folds(lm, avg_lm, 'AAC')
p2 <- plot_folds(lasso, avg_lasso, 'AAC')
p3 <- plot_folds(ridge, avg_ridge, 'AAC')
p4 <- plot_folds(en, avg_en, 'AAC')
p5 <- plot_folds(rf, avg_rf, 'AAC')

png("SupervisedSignatures/results/figures/AAC_indiv_folds.png", width=275, height=175, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, p5, nrow = 5, common.legend = TRUE, legend = "right")
dev.off()

### ===== Binarized Drug Response Model Performance Evaluation across Folds ===== ###

bin_path <- "SupervisedSignatures/results/data/bin/"

# load in indiv pset results
log <- read.csv(paste0(bin_path, "log/log.csv"))
lasso <- read.csv(paste0(bin_path, "lasso/lasso.csv"))
ridge <- read.csv(paste0(bin_path, "ridge/ridge.csv"))
en <- read.csv(paste0(bin_path, "en/en.csv"))
rf <- read.csv(paste0(bin_path, "rf/rf.csv"))       

# load in averaged results
avg_log <- read.csv(paste0(bin_path, "averaged/log.csv"))
avg_lasso <- read.csv(paste0(bin_path, "averaged/lasso.csv"))
avg_ridge <- read.csv(paste0(bin_path, "averaged/ridge.csv"))
avg_en <- read.csv(paste0(bin_path, "averaged/en.csv"))
avg_rf <- read.csv(paste0(bin_path, "averaged/rf.csv"))

# plot results across folds
p1 <- plot_folds(log, avg_log, 'bin')
p2 <- plot_folds(lasso, avg_lasso, 'bin')
p3 <- plot_folds(ridge, avg_ridge, 'bin')
p4 <- plot_folds(en, avg_en, 'bin')
p5 <- plot_folds(rf, avg_rf, 'bin')

png("SupervisedSignatures/results/figures/bin_indiv_folds.png", width=275, height=175, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, p5, nrow = 5, common.legend = TRUE, legend = "right")
dev.off()