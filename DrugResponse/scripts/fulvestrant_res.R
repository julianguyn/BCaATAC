# Plot scatter plot of Fulvestrant response + alternative drugs

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
})

###########################################################
# Load in data
###########################################################

# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")

# load in drug sensitivity data
load("DrugResponse/results/data/sensitivity_data.RData")
gray_sen <- gray_sen[,which(colnames(gray_sen) %in% samples$sample)] 
ctrp_sen <- ctrp_sen[,which(colnames(ctrp_sen) %in% samples$sample)] 
ccle_sen <- ccle_sen[,which(colnames(ccle_sen) %in% samples$sample)]
gdsc_sen <- gdsc_sen[,which(colnames(gdsc_sen) %in% samples$sample)]

# read in signature scores
signature_scores <- read.table("Signatures/data/bca_sign.Zscore.txt")
# remove duplicates
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

signature_scores <- signature_scores[,which(colnames(signature_scores) %in% samples$file)]
colnames(signature_scores) <- samples[match(colnames(signature_scores), samples$file), ]$sample
rownames(signature_scores) <- paste0("Signature", 1:6)

###########################################################
# List biomarkers of interest
###########################################################

biomarkers <- data.frame(signature = rep("Signature4", 2),
                         drug = c("Fulvestrant", "Teniposide"))

###########################################################
# Keep only drugs of interest
###########################################################

# function to subset for drugs
subset_drugs <- function(drug_sen, label) {
    drugs <- biomarkers$drug
    drug_sen <- drug_sen[rownames(drug_sen) %in% drugs,]
    if (nrow(drug_sen) > 0) {
        drug_sen$drug <- rownames(drug_sen)
        res <- melt(drug_sen, id.vars = "drug")
        colnames(res) <- c("Drug", "Sample", "AAC")
        res$PSet <- label
        return(res)
    }
    return(data.frame(matrix(nrow=0, ncol=2)))
}

ubr1_sen <- subset_drugs(ubr1_sen, "UHNBreast1")
ubr2_sen <- subset_drugs(ubr2_sen, "UHNBreast2")
gray_sen <- subset_drugs(gray_sen, "GRAY")
gcsi_sen <- subset_drugs(gcsi_sen, "gCSI")
gdsc_sen <- subset_drugs(gdsc_sen, "GDSC")
ccle_sen <- subset_drugs(ccle_sen, "CCLE")
ctrp_sen <- subset_drugs(ctrp_sen, "CTRP")

###########################################################
# Create dataframe for plotting
###########################################################

# merge drug response data
toPlot <- rbind(ubr1_sen, ubr2_sen, gray_sen, gcsi_sen, 
                gdsc_sen, ccle_sen, ctrp_sen)

# get signature scores
toPlot$Score <- NA
for (j in 1:nrow(toPlot)) { 
    score = signature_scores[rownames(signature_scores) == "Signature4",colnames(signature_scores) == toPlot$Sample[j]]
    if (is.numeric(score)) {
        toPlot$Score[j] <-  score
    } else {
        toPlot$Score[j] <- NA
    }
}

# compute correlation between siganture score and aac for each pset
corr_res <- data.frame(matrix(nrow = 0, ncol = 2))
toPlot$pset_drug <- paste(toPlot$PSet, toPlot$Drug, sep = "-")
for (pset_drug in unique(toPlot$pset_drug)) {
    tmp <- toPlot[toPlot$pset_drug == pset_drug,]
    corr <- cor(tmp$AAC, tmp$Score, use="complete.obs", method = "pearson")

    corr_res <- rbind(corr_res, data.frame(PSet_Drug = pset_drug, Correlation = corr))
}

###########################################################
# Plot Scatterplot
###########################################################

# set up palette for plotting
pal <- c("CTRP-Fulvestrant" = "#7D2E68", "GDSC-Fulvestrant" = "#C97B84", 
         "CTRP-Teniposide" = "#4D6D9D", "GDSC-Teniposide" = "#97C0CE")

# set pset_drug factor
toPlot$pset_drug <- factor(toPlot$pset_drug, 
    levels = c("CTRP-Fulvestrant", "GDSC-Fulvestrant", 
               "CTRP-Teniposide", "GDSC-Teniposide"))
   
# get x axis limits for plotting
x <- max(max(toPlot$Score), abs(min(toPlot$Score)))

# scatter plot
png("DrugResponse/results/figures/fulR/fulR.png", width=160, height=125, units='mm', res = 600, pointsize=80)
ggplot(toPlot, aes(x = Score, y = AAC, fill = pset_drug)) + geom_point(size = 3, shape = 21) + 
    geom_smooth(method = "lm", se=F, aes(color = pset_drug)) + 
    scale_fill_manual(values = pal[names(pal) %in% unique(toPlot$pset_drug)]) + 
    scale_color_manual(values = pal[names(pal) %in% unique(toPlot$pset_drug)]) +
    guides(color = 'none', fill = guide_legend(override.aes=list(values = pal[names(pal) %in% unique(toPlot$pset_drug)], linetype = 0))) +
    theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                            legend.key.size = unit(0.7, 'cm')) +
    labs(x = "Signature 4 Similarity Score", fill = "PSet-Drug Pair")
dev.off()
