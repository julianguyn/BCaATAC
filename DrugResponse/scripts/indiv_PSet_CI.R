setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(survcomp)
    library(wesanderson)
    library(ggplot2)
    library(ggh4x)
    library(reshape2)
})

###########################################################
# Load in data
###########################################################

# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")

# load in signature
signature_scores <- read.table("Signatures/data/bca_sign.Zscore.txt")
# remove duplicates
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

signature_scores <- signature_scores[,which(colnames(signature_scores) %in% samples$file)]


###########################################################
# Initiate dataframe to store signature associations
###########################################################

# data frame to keep track of number of significant associations
sig_com <- as.data.frame(matrix(ncol = 7, nrow = 0))
colnames(sig_com) <- c("pset", "signature", "drug", "ci", "pvalue", "FDR", "label")


###########################################################
# Define functions to compute concordance index
###########################################################

# function to compute CI
computeCI <- function(signature_scores, sensitivity_data, label) {

    # create data frame to hold results
    combinations <- as.data.frame(matrix(data = NA, nrow = nrow(signature_scores) * nrow(sensitivity_data), ncol = 7))
    colnames(combinations) <- c("signature", "drug", "ci", "pvalue", "se", "upper", "lower")
    combinations$signature <- rep(rownames(signature_scores), nrow(sensitivity_data))
    combinations$drug <- rep(rownames(sensitivity_data), each = nrow(signature_scores))

    # compute concordance index
    for (i in 1:nrow(combinations)){
        print(paste0(i, " out of ", nrow(combinations), " complete"))

        ci <- survcomp::concordance.index(as.numeric(sensitivity_data[combinations[,2][i],]), 
                                            # row from sensitivity_data with the drug for all samples
                                            surv.time = as.numeric(unlist(-signature_scores[combinations[,1][i],])), 
                                            # row of signature i in signature scores matrix with cell lines as columns
                                            surv.event = rep(1,length(sensitivity_data)), 
                                            # df of all drugs as rows with all samples as columns
                                            outx = TRUE, method="noether", na.rm = TRUE)

        combinations$pvalue[i] <- ci$p.value
        combinations$ci[i] <- ci$c.index
        combinations$se[i] <- ci$se
        combinations$upper[i] <- ci$upper
        combinations$lower[i] <- ci$lower
    }

    # filtering and multiple test correction
    combinations <- combinations[complete.cases(combinations$pvalue),]
    combinations$FDR <- p.adjust(combinations$pvalue, method = "BH", n = length(combinations$pvalue))
    combinations$FDRsig <- ifelse(combinations$FDR < 0.05, TRUE, FALSE)

    # rename signatures
    combinations$signature <- paste0("Signature", as.numeric(gsub("sign", "", combinations$signature)) + 1)

    # format dataframe for plotting (order by ci and add variable rank)
    combinations <- combinations[order(combinations$ci),]
    combinations$rank <- 1:nrow(combinations)
    combinations$pairs <- paste0(combinations$signature, "_", combinations$drug)
    combinations$pset <- c(rep(label, nrow(combinations)))
    
    return(combinations)
}

# function to save significant associations
saveSig <- function(sig_com, combinations, label) {

    # save signatures associated with drug sensitivity:
    sen <- combinations[which(combinations$ci > 0.6 & combinations$FDR < 0.05),]
    sig_com <- rbind(sig_com, data.frame(pset = c(rep(label, nrow(sen))), signature = sen$signature, drug = sen$drug, ci = sen$ci, pvalue = sen$pvalue, FDR = sen$FDR))
    write.csv(sen[order(sen$FDR, -sen$ci),], file = paste0("DrugResponse/results/tables/indiv_PSet_CI/", label, "_sen.csv"), row.names = F)

    # save signatures associated with drug resistance:
    res <- combinations[which(combinations$ci < 0.4 & combinations$FDR < 0.05),]
    sig_com <- rbind(sig_com, data.frame(pset = c(rep(label, nrow(res))), signature = res$signature, drug = res$drug, ci = res$ci, pvalue = res$pvalue, FDR = res$FDR))
    write.csv(res[order(res$FDR, res$ci),], file = paste0("DrugResponse/results/tables/indiv_PSet_CI/", label, "_res.csv"), row.names = F)

    return(sig_com)
}


###########################################################
# Load in PSet sensitivity
###########################################################

# UHN Breast PSet
uhnbreast <- readRDS("DrugResponse/data/PSet_UHNBreast.rds") |> updateObject()
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(uhnbreast, sensitivity.measure = "aac_recomputed",  fill.missing = F))
ubr1_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #43 unique cell lines remain, 7 dups
rm(uhnbreast)

# UHN Breast2 PSet
uhnbreast2 <- readRDS("DrugResponse/data/PharmacoSet.RDS")
sensitivity_data <- dcast(uhnbreast2@treatmentResponse$profiles, treatmentid ~ sampleid, value.var = "aac_recomputed")
rownames(sensitivity_data) <- sensitivity_data$treatmentid
sensitivity_data <- sensitivity_data[,-c(1)]
# mapping cell line names
mapping <- c("BT20" = "BT-20", "BT474" = "BT-474", "BT549" = "BT-549", "CAL120" = "CAL-120", "CAL148" = "CAL-148", 
             "CAL51" = "CAL-51", "CAMA1" = "CAMA-1", "EFM19" = "EFM-19", "HDQP1" = "HDQ-P1", "HS578T" = "Hs 578T",
             "JIMT1" = "JIMT-1", "KPL1" = "KPL-1", "LY2" = "MCF-7/LY2", "MCF7" = "MCF-7", "MDAMB157" = "MDA-MB-157",
             "MDAMB175VII" = "MDA-MB-175-VII", "MDAMB231" = "MDA-MB-231", "MDAMB361" = "MDA-MB-361", 
             "MDAMB436" = "MDA-MB-436", "MDAMB468" = "MDA-MB-468", "MX1" = "MX-1", "SKBR3" = "SK-BR-3", 
             "SKBR5" = "SK-BR-5", "SKBR7" = "SK-BR-7", "SUM149" = "SUM149PT", "SUM159" = "SUM159PT", "T47D" = "T-47D",
             "UACC812" = "UACC-812", "ZR751" = "ZR-75-1")
for (i in 1:length(colnames(sensitivity_data))) {
    cell = colnames(sensitivity_data)[i]
    if (cell %in% names(mapping)) {colnames(sensitivity_data)[i] <- unname(mapping[cell])}
}
# mapping drug names
mapping <- c("carboplatin" = "Carboplatinum", "dipyridamole" = "Dipyridamole", 
             "epirubicin" = "Epirubicin", "eribulin" = "Eribulin", "fluvastatin" = "Fluvastatin", 
             "herceptin" = "Herceptin", "lapatinib" = "Lapatinib", "paclitaxel" = "Paclitaxel")
for (i in 1:length(rownames(sensitivity_data))) {
    cell = rownames(sensitivity_data)[i]
    if (cell %in% names(mapping)) {rownames(sensitivity_data)[i] <- unname(mapping[cell])}
}
ubr2_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #42 unique cell lines remain
rm(uhnbreast2)

# Gray PSet
gray <- readRDS("DrugResponse/data/PSet_GRAY2017.rds") |> updateObject()
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(gray, sensitivity.measure = "aac_recomputed",  fill.missing = F))
gray_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #42 unique cell lines remain
rm(gray)

# gCSI
gcsi <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/gCSI.rds") |> updateObject()
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(gcsi, sensitivity.measure = "aac_recomputed",  fill.missing = F))
gcsi_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #25 unique cell lines remain, 7 dups
rm(gcsi)

# GDSC
gdsc <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/GDSC2-8.2.rds") |> updateObject()
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed",  fill.missing = F))
gdsc_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #34 unique cell lines remain, 7 dups
rm(gdsc)

# CTRP
load("C:/Users/julia/Desktop/ncRNA Code/data/PSets/CTRP.rds")
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(CTRP, sensitivity.measure = "aac_recomputed",  fill.missing = F))
ctrp_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #32 unique cell lines remain, 7 dups
rm(CTRP)

# CCLE
ccle <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/CCLE.rds") |> updateObject()
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "aac_recomputed",  fill.missing = F))
ccle_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #23 unique cell lines remain, 7 dups
rm(ccle)


###########################################################
# Subset signature associations
###########################################################

# UHN Breast
ubr1_sam <- samples[which(samples$sample %in% colnames(ubr1_sen)),]
ubr1_sig <- signature_scores[,which(colnames(signature_scores) %in% ubr1_sam$file)]

# UHN Breast2
ubr2_sam <- samples[which(samples$sample %in% colnames(ubr2_sen)),]
ubr2_sig <- signature_scores[,which(colnames(signature_scores) %in% ubr2_sam$file)]

# GRAY
gray_sam <- samples[which(samples$sample %in% colnames(gray_sen)),]
gray_sig <- signature_scores[,which(colnames(signature_scores) %in% gray_sam$file)]

# gCSI
gcsi_sam <- samples[which(samples$sample %in% colnames(gcsi_sen)),]
gcsi_sig <- signature_scores[,which(colnames(signature_scores) %in% gcsi_sam$file)]

# GDSC
gdsc_sam <- samples[which(samples$sample %in% colnames(gdsc_sen)),]
gdsc_sig <- signature_scores[,which(colnames(signature_scores) %in% gdsc_sam$file)]

# CTRP
ctrp_sam <- samples[which(samples$sample %in% colnames(ctrp_sen)),]
ctrp_sig <- signature_scores[,which(colnames(signature_scores) %in% ctrp_sam$file)]

# CCLE
ccle_sam <- samples[which(samples$sample %in% colnames(ccle_sen)),]
ccle_sig <- signature_scores[,which(colnames(signature_scores) %in% ccle_sam$file)]


###########################################################
# Compute concordance index CAS-drug associations
###########################################################

# UHN Breast
ubr1_com <- computeCI(ubr1_sig, ubr1_sen, "uhnbreast1")
sig_com <- saveSig(sig_com, ubr1_com, "uhnbreast1")

# UHN Breast2
ubr2_com <- computeCI(ubr2_sig, ubr2_sen, "uhnbreast2")
sig_com <- saveSig(sig_com, ubr2_com, "uhnbreast2")

# GRAY
gray_com <- computeCI(gray_sig, gray_sen, "gray")
sig_com <- saveSig(sig_com, gray_com, "gray")

# gCSI
gcsi_com <- computeCI(gcsi_sig, gcsi_sen, "gcsi")
sig_com <- saveSig(sig_com, gcsi_com, "gcsi")

# GDSC
gdsc_com <- computeCI(gdsc_sig, gdsc_sen, "gdsc")
sig_com <- saveSig(sig_com, gdsc_com, "gdsc")

# CTRP
ctrp_com <- computeCI(ctrp_sig, ctrp_sen, "ctrp")
sig_com <- saveSig(sig_com, ctrp_com, "ctrp")

# CCLE
ccle_com <- computeCI(ccle_sig, ccle_sen, "ccle")
sig_com <- saveSig(sig_com, ccle_com, "ccle")


###########################################################
# Rename signatures
###########################################################

# functon to rename signatures
rename_sig <- function(com_df) {
    com_df$signature <- gsub("Signature", "CAS-", com_df$signature)
    com_df$pairs <- gsub("Signature", "CAS-", com_df$pairs)
    return(com_df)
}

sig_com$pairs <- paste(sig_com$signature, sig_com$drug, sep = "_")
sig_com <- rename_sig(sig_com)
ubr1_com <- rename_sig(ubr1_com)
ubr2_com <- rename_sig(ubr2_com)
gray_com <- rename_sig(gray_com)
gcsi_com <- rename_sig(gcsi_com)
gdsc_com <- rename_sig(gdsc_com)
ctrp_com <- rename_sig(ctrp_com)
ccle_com <- rename_sig(ccle_com)

###########################################################
# Save signature associations
###########################################################

save(sig_com, ubr1_com, ubr2_com, gray_com, gcsi_com, gdsc_com, ctrp_com, ccle_com, 
     file = "DrugResponse/results/data/indiv_PSet_CI.RData")


###########################################################
# Define palettes for plotting
###########################################################

# signature palette
pal <- c("CAS-1" = "#046C9A", "CAS-2" = "#BBADB9", "CAS-3" = "#7294D4", 
        "CAS-4" = "#E8E1D9", "CAS-5" = "#AFC5D8", "CAS-6" = "#DF9C93")
# sensitive vs resistant palette
pal2 <- c("#899DA4", "#BC4749")

###########################################################
# Plot waterfall plots
###########################################################

# TODO: turn this into a function for each pset (if needed)
png("DrugResponse/results/figures/indiv_PSet_CI/uhnbreast.png", width = 8, height = 5, res = 600, units = "in")
ggplot(ubr1_com, aes(x = rank, y = ci - 0.5, fill = ifelse(ci > 0.65 | ci < 0.35, signature, "Below CI Threshold"))) + #shift ci by 0.5 so that the baseline becomes 0.5
    geom_bar(stat="identity") + scale_y_continuous(labels = function(y) y + 0.5) + theme_classic() +
    scale_fill_manual(values = pal) +
    labs(fill = "Signature", y = "Concordance Index (CI)", x = "Signature-Drug Pairs") + 
    geom_hline(yintercept = c(-0.15, 0.15), linetype = "dotted") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggtitle("UHN Breast") + theme(plot.title = element_text(hjust = 0.5)) 
dev.off()


###########################################################
# Identify Class A Biomarkers
###########################################################

# filter by CI thresholds
strict_sig_com <- sig_com[which(sig_com$ci > 0.7 | sig_com$ci < 0.3),]
strict_sig_com <- strict_sig_com[strict_sig_com$FDR < 0.01,]

# order by CI
strict_sig_com <- strict_sig_com[order(strict_sig_com$ci),]
strict_sig_com$rank <- 1:nrow(strict_sig_com)

# format factors
strict_sig_com$rank <- as.factor(strict_sig_com$rank)
strict_sig_com$drug <- as.factor(strict_sig_com$drug)

# specify sensitivity and resistance
strict_sig_com$type <- ifelse(strict_sig_com$ci > 0.5, "Sensitivity", "Resistance")
strict_sig_com$type <- factor(strict_sig_com$type, levels = c("Sensitivity", "Resistance"))


###########################################################
# Plot Class A biomarkers
###########################################################

# plot vertical
png("DrugResponse/results/figures/indiv_PSet_CI/ClassA-vertical.png", width = 7, height = 15, res = 600, units = "in")
ggplot(strict_sig_com, aes(x = ci - 0.5, y = rank)) +
    geom_col(aes(fill = signature), color = "black") + scale_x_continuous(limits = c(-0.5, 0.5), labels = function(x) x + 0.5) +
    scale_y_discrete(breaks = strict_sig_com$rank, labels = strict_sig_com$drug) +
    scale_fill_manual(values = pal) +
    theme_classic() + geom_vline(xintercept = 0) + 
    theme(legend.text = element_text(size = 13),
          legend.title = element_text(size = 16),
          axis.text.x = element_text(size = 13),  
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) +
    labs(y = "Drug", title = "", x = "Concordance Index (CI)", fill = "CAS") 
dev.off()

# format for same graph but horizontal
strict_sig_com$rank <- nrow(strict_sig_com):1
strict_sig_com$rank <- as.factor(strict_sig_com$rank)
strict_sig_com$drug <- as.factor(strict_sig_com$drug)
 
# plot horizontal
png("DrugResponse/results/figures/indiv_PSet_CI/ClassA-horizonal.png", width = 11, height = 5, res = 600, units = "in")
ggplot(strict_sig_com, aes(x = ci - 0.5, y = rank)) +
    geom_col(aes(fill = signature), color = "black") + scale_x_continuous(limits = c(-0.5, 0.5), labels = function(x) x + 0.5) +
    scale_y_discrete(breaks = strict_sig_com$rank, labels = strict_sig_com$drug) +
    scale_fill_manual(values = pal) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) + geom_vline(xintercept = 0) +
    labs(y = "Drug", title = "", x = "Concordance Index (CI)", fill = "CAS") + coord_flip()
dev.off()


###########################################################
# Save Class A Biomarkers
###########################################################

write.csv(strict_sig_com, file = "DrugResponse/results/data/ClassA_Biomarkers.csv", 
          quote = F, row.names = F)

###########################################################
# Waterfall plot of Class A Biomarkers
###########################################################

# function to make plot for each signature
waterfallplot_signature <- function(signature) {
    # subset dataframe to keep only signature of interest
    toPlot <- strict_sig_com[strict_sig_com$signature == signature,]

    # order by CI
    toPlot <- toPlot[order(toPlot$ci),]
    toPlot$rank <- 1:nrow(toPlot)

    toPlot$rank <- as.factor(toPlot$rank)
    toPlot$drug <- as.factor(toPlot$drug)

    p <- ggplot(toPlot, aes(x = ci - 0.5, y = rank)) +
    geom_col(aes(fill = type), color = "black") + scale_x_continuous(limits = c(-0.5, 0.5), labels = function(x) x + 0.5) +
    scale_y_discrete(breaks = toPlot$rank, labels = toPlot$drug) + geom_vline(xintercept = 0) +
    scale_fill_manual("Association", values = pal2) + theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    labs(y = "Drug", title = signature, x = "Concordance Index (CI)") 

    return(p)
}

png("DrugResponse/results/figures/indiv_PSet_CI/signature/sig1.png", width = 5, height = 5, res = 600, units = "in")
waterfallplot_signature("Signature1")
dev.off()

png("DrugResponse/results/figures/indiv_PSet_CI/signature/sig2.png", width = 5, height = 5, res = 600, units = "in")
waterfallplot_signature("Signature2")
dev.off()

png("DrugResponse/results/figures/indiv_PSet_CI/signature/sig3.png", width = 5, height = 5, res = 600, units = "in")
waterfallplot_signature("Signature3")
dev.off()

png("DrugResponse/results/figures/indiv_PSet_CI/signature/sig4.png", width = 5, height = 5, res = 600, units = "in")
waterfallplot_signature("Signature4")
dev.off()

png("DrugResponse/results/figures/indiv_PSet_CI/signature/sig5.png", width = 5, height = 5, res = 600, units = "in")
waterfallplot_signature("Signature5")
dev.off()

png("DrugResponse/results/figures/indiv_PSet_CI/signature/sig6.png", width = 5, height = 5, res = 600, units = "in")
waterfallplot_signature("Signature6")
dev.off()


###########################################################
# Additional Class A plots
###########################################################

# count per signature
png("DrugResponse/results/figures/indiv_PSet_CI/count_per_sig.png", width = 6, height = 4, res = 600, units = "in")
ggplot(strict_sig_com, aes(x = signature, fill = type)) + 
  geom_bar(stat = "count", size = 0.5, position = "dodge", color = "black", width=0.6) +
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0)) +
  scale_fill_manual("Association", values = pal2) +
  theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  labs(x = "CAS", y = "Number of Associations")
dev.off()

# table of counts per signature
#table(strict_sig_com$signature, strict_sig_com$type)

# count per PSet
png("DrugResponse/results/figures/indiv_PSet_CI/count_per_pset.png", width = 9, height = 2, res = 600, units = "in")
ggplot(strict_sig_com, aes(x = "", fill = pset)) + 
  geom_bar(stat = "count", size = 0.5, position = "fill", color = "black", width=0.6) +
  scale_fill_manual("PSet", values = wes_palette("Cavalcanti1")) +
  theme_void() + coord_polar("y", start=0) +
  labs(x = "CAS", y = "Proportion of Associations") +
  facet_wrap(~ signature, nrow = 1)
dev.off()