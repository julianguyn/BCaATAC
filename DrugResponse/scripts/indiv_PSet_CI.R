setwd("C:/Users/julia/Documents/BCaATAC")

suppressMessages(library(PharmacoGx))
suppressMessages(library(survcomp))
library(ggplot2)
library(wesanderson)

# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")

# load in palette for plotting
pal <- wes_palette("Cavalcanti1")
pal2 <- wes_palette("GrandBudapest2")

# get signature scores
signature_scores <- read.table("Signatures/bca_sign.Zscore.txt")
signature_scores <- signature_scores[,which(colnames(signature_scores) %in% samples$file)]
# remove duplicates
#samples[samples$sample %in% samples$sample[duplicated(samples$sample)],]
remove <- c("HCC1395.", "BT549.", "MDA.MB.436.", "HCC1143.", "HCC3153.", "HCC1806.", "BT20.", "HCC1937.", "SUM149PT.", "SUM159PT.", "MDA.MB.231.", "HS578.", "HCC70.", "MDA.MB.468.")
samples <- samples[-which(samples$file %in% remove),]
signature_scores <- signature_scores[,which(colnames(signature_scores) %in% samples$file)]

# data frame to keep track of number of significant associations
sig_com <- as.data.frame(matrix(ncol = 6, nrow = 0))
colnames(sig_com) <- c("pset", "signature", "drug", "ci", "pvalue", "FDR")


### Function to compute CI ###
computeCI <- function(signature_scores, sensitivity_data) {

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

    # format dataframe for plotting (order by ci and add variable rank)
    combinations$signature <- gsub("sign", "signature", combinations$signature)
    combinations <- combinations[order(combinations$ci),]
    combinations$rank <- 1:nrow(combinations)
    combinations$pairs <- paste(combinations$signature, combinations$drug)

    return(combinations)
}

saveSig <- function(sig_com, combinations, label) {

    # save signatures associated with drug sensitivity:
    sen <- combinations[which(combinations$ci > 0.65),]
    sig_com <- rbind(sig_com, data.frame(pset = c(rep(label, nrow(sen))), signature = sen$signature, drug = sen$drug, ci = sen$ci, pvalue = sen$pvalue, FDR = sen$FDR))
    write.csv(sen[order(sen$FDR, -sen$ci),], file = paste0("DrugResponse/results/tables/indiv_PSet_CI/", label, "_sen.csv"), row.names = F)

    # save signatures associated with drug resistance:
    res <- combinations[which(combinations$ci < 0.35),]
    sig_com <- rbind(sig_com, data.frame(pset = c(rep(label, nrow(res))), signature = res$signature, drug = res$drug, ci = res$ci, pvalue = res$pvalue, FDR = res$FDR))
    write.csv(res[order(res$FDR, res$ci),], file = paste0("DrugResponse/results/tables/indiv_PSet_CI/", label, "_res.csv"), row.names = F)

    return(sig_com)
}

#######################
### UHN Breast PSet ###
#######################

# load in UHN Breast PSet
uhnbreast <- readRDS("DrugResponse/data/PSet_UHNBreast.rds")
uhnbreast <- updateObject(uhnbreast)

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(uhnbreast, sensitivity.measure = "aac_recomputed",  fill.missing = F))
ubr1_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #43 unique cell lines remain, 7 dups
rm(uhnbreast)

# subset signature scores to samples in dataset
ubr1_sam <- samples[which(samples$sample %in% colnames(ubr1_sen)),]
ubr1_sig <- signature_scores[,which(colnames(signature_scores) %in% ubr1_sam$file)]

# compute associations
ubr1_com <- computeCI(ubr1_sig, ubr1_sen)
sig_com <- saveSig(sig_com, ubr1_com, "uhnbreast")

# plot associations
png("DrugResponse/results/figures/indiv_PSet_CI/uhnbreast.png", width = 8, height = 5, res = 600, units = "in")
ggplot(ubr1_com, aes(x = rank, y = ci - 0.5, fill = ifelse(ci > 0.65 | ci < 0.35, signature, "Below CI Threshold"))) + #shift ci by 0.5 so that the baseline becomes 0.5
    geom_bar(stat="identity") + scale_y_continuous(labels = function(y) y + 0.5) + theme_classic() +
    scale_fill_manual(values = c("Other" = "grey", "signature2" = pal[1], "signature3" = pal[2], "signature4" = pal[5], "signature5" = pal[4])) +
    labs(fill = "Signature", y = "Concordance Index (CI)", x = "Signature-Drug Pairs") + 
    geom_hline(yintercept = c(-0.15, 0.15), linetype = "dotted") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggtitle("UHN Breast") + theme(plot.title = element_text(hjust = 0.5)) +
    #annotate(geom = "text", label = c("CI > 0.65", "CI < 0.35"), x = c(21.5, 21.5), y = c(0.18, -0.16), vjust = 1) +
    annotate(geom = "text", label = c("*", "*", "*"), x = c(1, 3, 43), y = c(-0.22, -0.17, 0.22), vjust = 1)
dev.off()


#################
### Gray PSet ###
#################

# load in Gray PSet
gray <- readRDS("DrugResponse/data/PSet_GRAY2017.rds")
gray <- updateObject(gray)

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(gray, sensitivity.measure = "aac_recomputed",  fill.missing = F))
gray_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #42 unique cell lines remain
rm(gray)

# subset signature scores to samples in dataset
gray_sam <- samples[which(samples$sample %in% colnames(gray_sen)),]
gray_sig <- signature_scores[,which(colnames(signature_scores) %in% gray_sam$file)]

# compute associations
gray_com <- computeCI(gray_sig, gray_sen)
sig_com <- saveSig(sig_com, gray_com, "gray")

# plot associations
png("DrugResponse/results/figures/indiv_PSet_CI/gray.png", width = 8, height = 5, res = 600, units = "in")
ggplot(gray_com, aes(x = rank, y = ci - 0.5, fill = ifelse(ci > 0.65 | ci < 0.35, signature, "Below CI Threshold"))) + #shift ci by 0.5 so that the baseline becomes 0.5
    geom_bar(stat="identity") + scale_y_continuous(labels = function(y) y + 0.5) + theme_classic() +
    scale_fill_manual(values = c("Other" = "grey", "signature0" = pal2[3], "signature1" = pal2[4], "signature2" = pal[1], "signature3" = pal[2], "signature4" = pal[5], "signature5" = pal[4])) +
    labs(fill = "Signature", y = "Concordance Index (CI)", x = "Signature-Drug Pairs") + 
    geom_hline(yintercept = c(-0.15, 0.15), linetype = "dotted") +
    ggtitle("GRAY") + theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()


#################
### gCSI PSet ###
#################

# load in gCSI PSet
gcsi <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/gCSI.rds")
gcsi <- updateObject(gcsi)

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(gcsi, sensitivity.measure = "aac_recomputed",  fill.missing = F))
gcsi_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #25 unique cell lines remain, 7 dups
rm(gcsi)

# subset signature scores to samples in dataset
gcsi_sam <- samples[which(samples$sample %in% colnames(gcsi_sen)),]
gcsi_sig <- signature_scores[,which(colnames(signature_scores) %in% gcsi_sam$file)]

# compute associations
gcsi_com <- computeCI(gcsi_sig, gcsi_sen)
sig_com <- saveSig(sig_com, gcsi_com, "gcsi")

# plot associations
png("DrugResponse/results/figures/indiv_PSet_CI/gcsi.png", width = 8, height = 5, res = 600, units = "in")
ggplot(gcsi_com, aes(x = rank, y = ci - 0.5, fill = ifelse(ci > 0.65 | ci < 0.35, signature, "Below CI Threshold"))) + #shift ci by 0.5 so that the baseline becomes 0.5
    geom_bar(stat="identity") + scale_y_continuous(labels = function(y) y + 0.5) + theme_classic() +
    scale_fill_manual(values = c("Other" = "grey", "signature1" = pal2[4], "signature2" = pal[1], "signature4" = pal[5])) +
    labs(fill = "Signature", y = "Concordance Index (CI)", x = "Signature-Drug Pairs") + 
    geom_hline(yintercept = c(-0.15, 0.15), linetype = "dotted") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggtitle("gCSI") + theme(plot.title = element_text(hjust = 0.5)) +
    annotate(geom = "text", label = c("*"), x = c(1), y = c(-0.18), vjust = 1)
dev.off()


################# 
### GDSC PSet ###
#################

# load in GDSC PSet
gdsc <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/GDSC2-8.2.rds")
gdsc <- updateObject(gdsc)

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed",  fill.missing = F))
gdsc_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #34 unique cell lines remain, 7 dups
rm(gdsc)

# subset signature scores to samples in dataset
gdsc_sam <- samples[which(samples$sample %in% colnames(gdsc_sen)),]
gdsc_sig <- signature_scores[,which(colnames(signature_scores) %in% gdsc_sam$file)]

# compute associations
gdsc_com <- computeCI(gdsc_sig, gdsc_sen)
sig_com <- saveSig(sig_com, gdsc_com, "gdsc")

# plot associations
png("DrugResponse/results/figures/indiv_PSet_CI/gdsc.png", width = 8, height = 5, res = 600, units = "in")
ggplot(gdsc_com, aes(x = rank, y = ci - 0.5, fill = ifelse(ci > 0.65 | ci < 0.35, signature, "Below CI Threshold"))) + #shift ci by 0.5 so that the baseline becomes 0.5
    geom_bar(stat="identity") + scale_y_continuous(labels = function(y) y + 0.5) + theme_classic() +
    scale_fill_manual(values = c("Other" = "grey", "signature0" = pal2[3], "signature1" = pal2[4], "signature2" = pal[1], "signature3" = pal[2], "signature4" = pal[5], "signature5" = pal[4])) +
    labs(fill = "Signature", y = "Concordance Index (CI)", x = "Signature-Drug Pairs") + 
    geom_hline(yintercept = c(-0.15, 0.15), linetype = "dotted") +
    ggtitle("GDSC2") + theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
dev.off()


#################
### CTRP PSet ###
#################

# load in CTRP PSet
load("C:/Users/julia/Desktop/ncRNA Code/data/PSets/CTRP.rds")

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(CTRP, sensitivity.measure = "aac_recomputed",  fill.missing = F))
ctrp_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #32 unique cell lines remain, 7 dups
rm(CTRP)

# subset signature scores to samples in dataset
ctrp_sam <- samples[which(samples$sample %in% colnames(ctrp_sen)),]
ctrp_sig <- signature_scores[,which(colnames(signature_scores) %in% ctrp_sam$file)]

# compute associations
ctrp_com <- computeCI(ctrp_sig, ctrp_sen)
sig_com <- saveSig(sig_com, ctrp_com, "ctrp")

# plot associations
png("DrugResponse/results/figures/indiv_PSet_CI/ctrp.png", width = 8, height = 5, res = 600, units = "in")
ggplot(ctrp_com, aes(x = rank, y = ci - 0.5, fill = ifelse(ci > 0.65 | ci < 0.35, signature, "Below CI Threshold"))) + #shift ci by 0.5 so that the baseline becomes 0.5
    geom_bar(stat="identity") + scale_y_continuous(labels = function(y) y + 0.5) + theme_classic() +
    scale_fill_manual(values = c("Other" = "grey", "signature0" = pal2[3], "signature1" = pal2[4], "signature2" = pal[1], "signature3" = pal[2], "signature4" = pal[5], "signature5" = pal[4])) +
    labs(fill = "Signature", y = "Concordance Index (CI)", x = "Signature-Drug Pairs") + 
    geom_hline(yintercept = c(-0.15, 0.15), linetype = "dotted") +
    ggtitle("CTRP") + theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
dev.off()


#################
### CCLE PSet ###
#################

# load in CCLE PSet
ccle <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/CCLE.rds")
ccle <- updateObject(ccle)

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "aac_recomputed",  fill.missing = F))
ccle_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #23 unique cell lines remain, 7 dups
rm(ccle)

# subset signature scores to samples in dataset
ccle_sam <- samples[which(samples$sample %in% colnames(ccle_sen)),]
ccle_sig <- signature_scores[,which(colnames(signature_scores) %in% ccle_sam$file)]

# compute associations
ccle_com <- computeCI(ccle_sig, ccle_sen)
sig_com <- saveSig(sig_com, ccle_com, "ccle")

# plot associations
png("DrugResponse/results/figures/indiv_PSet_CI/ccle.png", width = 8, height = 5, res = 600, units = "in")
ggplot(ccle_com, aes(x = rank, y = ci - 0.5, fill = ifelse(ci > 0.65 | ci < 0.35, signature, "Below CI Threshold"))) + #shift ci by 0.5 so that the baseline becomes 0.5
    geom_bar(stat="identity") + scale_y_continuous(labels = function(y) y + 0.5) + theme_classic() +
    scale_fill_manual(values = c("Other" = "grey", "signature0" = pal2[3], "signature1" = pal2[4], "signature2" = pal[1], "signature4" = pal[5])) +
    labs(fill = "Signature", y = "Concordance Index (CI)", x = "Signature-Drug Pairs") + 
    geom_hline(yintercept = c(-0.15, 0.15), linetype = "dotted") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggtitle("CCLE") + theme(plot.title = element_text(hjust = 0.5)) +
    annotate(geom = "text", label = c("*", "*", "*", "*"), x = c(1, 3.2, 141.8, 143), y = c(-0.208, -0.175, 0.209, 0.222), vjust = 1) 
dev.off()



####################################
### All Significant Associations ###
####################################

### Count per Signature ###
sig_com$type <- ifelse(sig_com$ci > 0.5, "Sensitivity", "Resistance")

png("DrugResponse/results/figures/indiv_PSet_CI/count_per_sig.png", width = 6, height = 4, res = 600, units = "in")
ggplot(sig_com, aes(x = signature, fill = type)) + 
  geom_bar(stat = "count", size = 0.5, position = "dodge", color = "black", width=0.6) +
  scale_y_continuous(limits = c(0, 75), expand = c(0, 0)) +
  scale_fill_manual("Association", values = c("Sensitivity" = pal[1], "Resistance" = pal[5])) +
  theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  labs(x = "Signature", y = "Number of Associations")
dev.off()


### Count per Drug ###
drug_count <- as.data.frame(table(sig_com$drug)[order(-table(sig_com$drug))])
colnames(drug_count) <- c("Drug", "Count")

# keep only drugs with count > 2
drug_count <- drug_count[drug_count$Count > 2,]

# get number of psets drug counts span
drug_count$PSet = 0
for (i in 1:length(drug_count$Drug)) {
    
    # subset all significant combinations to those with the drug
    subsetdf <- sig_com[sig_com$drug == drug_count$Drug[i],]

    # get number of unique psets
    pset_count <- length(unique(subsetdf$pset))
    drug_count$PSet[i] <- pset_count
}
drug_count$PSet <- as.factor(drug_count$PSet)

png("DrugResponse/results/figures/indiv_PSet_CI/count_per_drug.png", width = 6, height = 12, res = 600, units = "in")
ggplot(drug_count, aes(x = reorder(Drug, -as.numeric(factor(Drug))), y = Count, fill = PSet)) + 
  geom_col(size = 0.5, position = "identity", color = "black") +
  coord_flip() + theme_classic() +
  scale_fill_manual("Number of PSets", values = c("1" = pal2[2], "2" = pal2[4], "4" = pal2[3])) +
  labs(x = "Drug", y = "Count")
dev.off()


### Signature-Drug pairs significant in >1 PSet ###
sig_com$pairs <- paste0(sig_com$signature, sig_com$drug)
multiple <- sig_com$pairs[duplicated(sig_com$pairs)]

df <- sig_com[which(sig_com$pairs %in% multiple),]
df$pset <- toupper(df$pset)

png("DrugResponse/results/figures/indiv_PSet_CI/double_pairs.png", width = 10, height = 5, res = 600, units = "in")
ggplot(df, aes(x = pset, y = ci, fill = signature)) + 
  geom_bar(stat = "identity", size = 0.5, position = position_dodge(width = 0.9), color = "black") + 
  facet_grid(~ as.factor(drug), scales = "free_x") +
  scale_y_continuous(limits = c(0, 0.8), expand=c(0,0)) +
  scale_fill_manual("Signature", values = c("signature0" = pal2[3], "signature1" = pal2[4], "signature3" = pal[2], "signature4" = pal[5])) +
  theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  labs(x = "PSet", y = "Concordance Index (CI)")
dev.off()

save(sig_com, ubr1_com, gray_com, gcsi_com, gdsc_com, ctrp_com, ccle_com, file = "DrugResponse/results/data/indiv_PSet_CI.RData")