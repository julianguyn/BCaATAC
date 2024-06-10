setwd("C:/Users/julia/Documents/BCaATAC")

library(readxl)
library(ggplot2)
library(RColorBrewer)

# read in drug response data
response <- as.data.frame(read_excel("DrugResponsePDX/data/drugresponse/DrugResponse_PDX.xlsx", sheet = 1))

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
strict_sig_com$drug <- as.factor(strict_sig_com$drug)

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


# find overlaps with Dave's group
grep("263", X)
grep("5305", X)
grep("8205", X)
grep("673", X)
grep("011", X)
grep("400945", X)
grep("402257", X)
grep("5461", X)
grep("G9A", X)
grep("400", X)
grep("0642", X)
grep("PQ", X)
grep("carbop", X)
grep("potomab", X)
grep("eribulin", X)
grep("everolimus", X)
grep("palbociclub", X)
grep("tax", X)
grep("zumab", X)
grep("fluv", X)