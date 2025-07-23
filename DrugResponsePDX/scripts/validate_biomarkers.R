setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(readxl)
    library(ggplot2)
    library(RColorBrewer)
    library(dplyr)
})

source("DrugResponsePDX/scripts/source/helper.R")

###########################################################
# Load in data
###########################################################

# read in drug response data
old_response <- as.data.frame(read_excel("DrugResponsePDX/data/drugresponse/DrugResponse_PDX.xlsx", sheet = 1))
new_xeva_auc <- read.csv("DrugResponsePDX/data/Jul232025-Xeva/auc_reps.csv")
new_xeva_mRECIST <- read.csv("DrugResponsePDX/data/Jul232025-Xeva/mRECIST_reps.csv")

# read in signature scores
scores <- as.data.frame(t(read.table("DrugResponsePDX/data/chromvar/bca_sign.Zscore.txt")))
colnames(scores) <- paste0("Signature", 1:6)


###########################################################
# Format sample names
###########################################################

# using sourced map_pdx() to standardize PDX names
rownames(scores) <- map_pdx(rownames(scores))
old_response$patient.id <- map_pdx(old_response$patient.id)
new_xeva_auc$patient.id <- map_pdx(new_xeva_auc$patient.id)
new_xeva_mRECIST$patient.id <- map_pdx(new_xeva_mRECIST$patient.id)


###########################################################
# Assign signature scores
###########################################################

# using sourced get_ARCHE() to extract PDX ARCHE scores
old_response <- get_ARCHE(old_response)
new_xeva_auc <- get_ARCHE(new_xeva_auc)
new_xeva_mRECIST <- get_ARCHE(new_xeva_mRECIST)


###########################################################
# Assess ARCHE treatment response associations 
###########################################################

# using sourced 
assess_ARCHE_mRECIST <- function(df, drug) {
    # subset for variables of interest
    df <- df[df$drug == drug,]
    if (c('slope') %in% colnames(df)) {
        df[,-which(colnames(df) %in% c('slope', 'AUC'))]
    }
    df <- df[complete.cases(df),]
    df$mRECIST <- factor(df$mRECIST, levels = c("CR", "PR", "SD", "PD"))

    # iterate through each ARCHE
    for (arche in paste0('ARCHE', 1:6)) {
        p1 <- waterfall_mRECIST(df, arche, drug)
    }
}



# df <- old_response[old_response$drug == "ERIBULIN",]


###########################################################
# Plot waterfall by AUC
###########################################################

# order samples
pac$AUC <- as.numeric(pac$AUC)
pac_auc <- pac[order(pac$AUC, decreasing = T),]
pac_auc$rank <- 1:nrow(pac_auc)

# plot waterfall plot coloured by signature score
png("DrugResponsePDX/results/figures/paclitaxel_pdx_sig_waterfall.png", width=175, height=125, units='mm', res = 600, pointsize=80)
ggplot(pac_auc, aes(x = rank, y = AUC, fill = sig5)) + 
    geom_bar(stat = "identity", color = "black") + geom_hline(yintercept = 0) +
    scale_fill_gradientn(colors = c("#BC4749", "#F8F1F8", "#077293"), limits = c(-45, 45)) +
    theme_classic() + theme(legend.key.size = unit(0.8, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = "PDX Model", y = "AUC", fill = "Signature 5\nScore") 
dev.off()


###########################################################
# Plot waterfall by slope
###########################################################

# order samples
pac$slope <- as.numeric(pac$slope)
pac_slp <- pac[order(pac$slope, decreasing = T),]
pac_slp$rank <- 1:nrow(pac_slp)

# plot waterfall plot coloured by signature score
png("DrugResponsePDX/results/figures/paclitaxel_pdx_slope_waterfall.png", width=175, height=125, units='mm', res = 600, pointsize=80)
ggplot(pac_slp, aes(x = rank, y = slope, fill = sig5)) + 
    geom_bar(stat = "identity", color = "black") + geom_hline(yintercept = 0) +
    scale_fill_gradientn(colors = c("#BC4749", "#F8F1F8", "#077293"), limits = c(-45, 45)) +
    theme_classic() + theme(legend.key.size = unit(0.8, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = "PDX Model", y = "Slope", fill = "Signature 5\nScore") 
dev.off()


###########################################################
# Plot waterfall of average signature score per mRECIST
###########################################################

# average signature 5 score for each mRECIST category
pac_mr <- pac_mr %>%
  group_by(mRECIST) %>%
  summarize(mean_score = mean(sig5), sd_score = sd(sig5))


png("DrugResponsePDX/results/figures/paclitaxel_pdx_mrecist.png", width = 4, height = 5, res = 600, units = "in")
ggplot(pac_mr, aes(x = mRECIST, y = mean_score)) +
  geom_bar(stat = "identity", color = "black", aes(fill = mRECIST)) + ylim(c(-50, 50)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score), width = 0.2) +
  scale_fill_manual(values = c("#136F63", "#9DCBBA", "#FFADA1", "#B02E0C"),
                      labels = c('CR' = 'Complete\nResponse', 
                              'PD' = 'Progressive\nDisease', 
                              'SD' = 'Stable\nDisease',
                              'PR' = 'Partial\nResponse')) +
  scale_x_discrete(labels = c('CR' = 'Complete\nResponse', 
                              'PD' = 'Progressive\nDisease', 
                              'SD' = 'Stable\nDisease',
                              'PR' = 'Partial\nResponse')) +
  labs(x = "\nmRECIST", y = "Signature 5 Similarity Score") +
  theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5), legend.position = "none")
dev.off()


###########################################################
# Plot waterfall of prediction accuracy
###########################################################

# keep only AUC
response <- response[,colnames(response) %in% c("patient.id", "AUC", "drug")]
#dcast(response, drug ~ patient.id, value.var = "AUC") cant' do this bc not all samples have all drugs


# read in signature scores
scores <- as.data.frame(t(read.table("DrugResponsePDX/data/chromvar/bca_sign.Zscore.txt")))
colnames(scores) <- paste0("Signature", 1:6)
# sample colname


# read in preclinical biomarkers
classA <- read.csv("DrugResponse/results/data/ClassA_Biomarkers.csv")
classB <- read.csv("DrugResponse/results/data/ClassB_Biomarkers.csv")
classC <- read.csv("DrugResponse/results/data/ClassC_Biomarkers.csv")

# get all drugs with predictions
all_drugs <- unique(c(classA$drug, classB$drug, classC$drug))

# TODO: remove "pitstop2"

# subset for just paclitaxel
pac <- response[response$drug == "TAXOL",]
pac <- pac[-which(is.na(pac$patient.id)),]
pac$AUC <- as.numeric(pac$AUC)
pac$sig5 <- NA

# get signature scores
for (i in 1:nrow(pac)) {
    sample = pac$patient.id[i]
    pac$sig5[i] <- scores[gsub("_", "", gsub("X", "", rownames(scores))) == sample,]$Signature5
}

pac$quadrant <- ifelse(pac$AUC > 0, ifelse(pac$sig5 > 0, "Q1", "Q2"), ifelse(pac$sig5 < 0.5, "Q3", "Q4"))

# plot 
x = max(abs(min(pac$sig5)), max(pac$sig5))
y = max(abs(min(pac$AUC)), max(pac$AUC))

pac <- pac[order(pac$AUC, decreasing = T),]
pac$rank <- 1:nrow(pac)
pac$rank <- as.factor(pac$rank)

pac$pred <- ifelse(pac$quadrant %in% c("Q2", "Q3"), "Accurate", "Inaccurate") 
pac$pred <- factor(pac$pred, levels = c("Accurate", "Inaccurate"))

png("DrugResponsePDX/results/figures/paclitaxel_pdx_waterfall.png", width = 8, height = 5, res = 600, units = "in")
ggplot(pac, aes(x = AUC, y = rank)) +
    geom_col(aes(fill = pred), color = "black") + scale_x_continuous(limits = c(-y, y), labels = function(x) x + 0.5) +
    #scale_y_discrete(breaks = strict_sig_com$rank, labels = strict_sig_com$drug) +
    scale_fill_manual(values = c("#899DA4", "#BC4749")) + 
    theme_classic() + theme(axis.text.x = element_blank()) + geom_vline(xintercept = 0) +
    labs(y = "PDX Model", title = "", x = "True AUC", fill = "Signature Prediction") + coord_flip()
dev.off()


###########################################################
# Plot accuracy quadrants
###########################################################

png("DrugResponsePDX/results/figures/paclitaxel_pdx.png", width=160, height=125, units='mm', res = 600, pointsize=80)
ggplot(pac, aes(x = sig5, y = AUC, fill = patient.id)) + 
    geom_rect(fill = "#EDFFEF", color = NA, xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    geom_rect(fill = "#EDFFEF", color = NA, xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf) + 
    geom_rect(fill = "#FFECEC", color = "black", linetype = "dotted", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0) + 
    geom_rect(fill = "#FFECEC", color = "black", linetype = "dotted", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf) + 
    geom_point(size = 4, shape = 21) + scale_fill_brewer(palette = "Set3") +
    theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                            plot.title = element_text(hjust = 0.5, size = 16), legend.key.size = unit(0.7, 'cm')) +
    xlim(-x, x) + ylim(-y, y) + labs(x = "Signature 5 Similarity Score", fill = "PDX Model", title = "Paclitaxel") 
dev.off()


# subset for just UNC0642
unc <- response[response$drug == "UNC0642",]
unc <- unc[-which(is.na(unc$patient.id)),]
unc$AUC <- as.numeric(unc$AUC)
unc$sig4 <- NA

# get signature scores
for (i in 1:nrow(unc)) { unc$sig4[i] <- scores[gsub("_", "", gsub("X", "", rownames(scores))) == unc$patient.id[i],]$Signature4 }

# plot 
x = max(abs(min(unc$sig4)), max(unc$sig4))
y = max(abs(min(unc$AUC)), max(unc$AUC))

png("DrugResponsePDX/results/figures/unc0642_pdx.png", width=160, height=125, units='mm', res = 600, pointsize=80)
ggplot(unc, aes(x = sig4, y = AUC, fill = patient.id)) + 
    geom_rect(fill = "#EDFFEF", color = NA, xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    geom_rect(fill = "#EDFFEF", color = NA, xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf) + 
    geom_rect(fill = "#FFECEC", color = "black", linetype = "dotted", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0) + 
    geom_rect(fill = "#FFECEC", color = "black", linetype = "dotted", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf) + 
    geom_point(size = 4, shape = 21) + scale_fill_brewer(palette = "Set3") +
    theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                            plot.title = element_text(hjust = 0.5, size = 16), legend.key.size = unit(0.7, 'cm')) +
    xlim(-x, x) + ylim(-y, y) + labs(x = "Signature 4 Similarity Score", fill = "PDX Model", title = "UNC0642") 
dev.off()