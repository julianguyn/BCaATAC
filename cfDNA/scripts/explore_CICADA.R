# plot results from Sasha (top 10k sites in each CALS)

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(data.table)
})

###########################################################
# Load in data
###########################################################

# get list of files to open
path <- "cfDNA/data/2025-09-09_julia_signature_results/"
#files <- list.files(path)
#files <- files[grep("CICADA", files)]

# load in metadata
meta <- read.csv(paste0(path, "sample_metadata.csv"))

# get control samples
control <- meta$sample_id[meta$time_id == "On treatment control"]

# filter by TF
meta <- meta[meta$metrics_tf > 10,]

# get baseline and progression samples
baseline <- meta$sample_id[meta$time_id == "Baseline"]
baseline <- baseline[-which(baseline == "CICADA-0031-01")]
progress <- meta$sample_id[meta$time_id == "Progression"]
progress <- progress[-which(progress == "CICADA-0011-02")]

# get files for each group
b_files <- paste0(path, baseline, "/", baseline, ".GC_corrected.coverage.tsv")
p_files <- paste0(path, progress, "/", progress, ".GC_corrected.coverage.tsv")
c_files <- paste0(path, control, "/", control, ".GC_corrected.coverage.tsv")

###########################################################
# Function to get CALS score
###########################################################

# function to get CALS score from Griffin output matrix
keep <- c("sig1_top10k", "sig2_top10k", "sig3_top10k", "sig4_top10k", "sig5_top10k", "sig6_top10k")

score_CALS <- function(file, res) {
    griffin <- fread(file)
    griffin <- griffin[griffin$site_name %in% keep,]
    griffin <- data.frame(CALS = paste0("CALS-", 1:6),
                          Score = 1-griffin$central_coverage)
    res <- rbind(res, griffin)
    return(res)
}

###########################################################
# Get CALS score
###########################################################

# create dataframe to store results
res <- data.frame(matrix(nrow=0, ncol=2)) 
colnames(res) <- c("CALS", "Score")

b_res <- p_res <- c_res <- res

# get baseline results (30)
for (file in b_files) {
    b_res <- score_CALS(file, b_res)
}
b_res$label = "Baseline"

# get progression results (14)
for (file in p_files) {
    p_res <- score_CALS(file, p_res)
}
p_res$label = "Progression"

# get control results (3)
for (file in c_files) {
    c_res <- score_CALS(file, c_res)
}
c_res$label = "Control"


###########################################################
# Plot CALS scores
###########################################################

# create datafrme to plot
toPlot <- rbind(b_res, p_res, c_res)

pal <- c("#D07C7B", "#8E2C2C", "#F5C6C6")

# ai

library(dplyr)
library(ggplot2)

# Step 1: Summarize the data
summary_df <- toPlot %>%
  group_by(CALS, label) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    se_score = sd(Score, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# factor
summary_df$label <- factor(summary_df$label, levels = c("Baseline", "Progression", "Control"))

# Step 2: Plot
png("cfDNA/results/figures/CALS_score_CICADA.png", width=150, height=100, units='mm', res = 600, pointsize=80)
ggplot(summary_df, aes(x = CALS, y = mean_score, fill = label)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") + 
  geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                width = 0.2, position = position_dodge(width = 0.9)) +
  theme_classic() + 
  geom_hline(yintercept = 0) + 
  scale_fill_manual("Sample Type", values = pal) +
  theme(legend.key.size = unit(0.7, 'cm')) + 
  labs(x = "CALS", y = "CALS Score")
dev.off()
