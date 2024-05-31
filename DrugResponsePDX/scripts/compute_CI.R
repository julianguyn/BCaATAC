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