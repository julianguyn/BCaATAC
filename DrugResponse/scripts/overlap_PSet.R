### This script checks or overlapping cell lines and drugs across PSets
### Also correlations between PSet drug response for overlapping cell lines and drugs

setwd("C:/Users/julia/Documents/BCaATAC")

library(PharmacoGx)
library(ggplot2)
library(wesanderson)

# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")

#######################
### UHN Breast PSet ###
#######################

# load in UHN Breast PSet
uhnbreast <- readRDS("DrugResponse/data/PSet_UHNBreast.rds")
uhnbreast <- updateObject(uhnbreast)

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(uhnbreast, sensitivity.measure = "aac_recomputed",  fill.missing = F))
sensitivity_data <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #43 unique cell lines remain, 7 dups

# save drug sensitivity, cell lines and drugs
uhnbreast_cl <- colnames(sensitivity_data)
uhnbreast_dr <- rownames(sensitivity_data)


########################
### UHN Breast PSet2 ###
########################

# load in UHN Breast PSet
load("DrugResponse/data/PharmacoGxCards/data/TFRI_TNBC_UHN_All_PSet.RData")

# save cell lines and drugs
# manually add into PSet
# DU-4475 = DU4475
# Hs-578-T = Hs 578T
# MCF10A = MCF-10A
# MCF7 = MCF-7
# LY2 = MCF-7/LY2
# SUM 149PT = SUM149PT
# SUM 159PT = SUM159PT
# SUM 52PE = SUM52PE
# T47D = T-47D
uhnbreast2_cl <- rownames(TFRI_TNBC_UHN@curation$cell)
uhnbreast2_cl <- unique(samples$sample)
uhnbreast2_dr <- rownames(TFRI_TNBC_UHN@curation$drug)

#################
### Gray PSet ###
#################

# load in Gray PSet
gray <- readRDS("DrugResponse/data/PSet_GRAY2017.rds")
gray <- updateObject(gray)

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(gray, sensitivity.measure = "aac_recomputed",  fill.missing = F))
sensitivity_data <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #42 unique cell lines remain


# save drug sensitivity, cell lines and drugs
gray_cl <- colnames(sensitivity_data)
gray_dr <- rownames(sensitivity_data)


#################
### gCSI PSet ###
#################

# load in gCSI PSet
gcsi <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/gCSI.rds")
gcsi <- updateObject(gcsi)

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(gcsi, sensitivity.measure = "aac_recomputed",  fill.missing = F))
sensitivity_data <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #25 unique cell lines remain, 7 dups

# save drug sensitivity, cell lines and drugs
gcsi_cl <- colnames(sensitivity_data)
gcsi_dr <- rownames(sensitivity_data)

#################
### GDSC PSet ###
#################

# load in GDSC PSet
gdsc <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/GDSC2-8.2.rds")
gdsc <- updateObject(gdsc)

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed",  fill.missing = F))
sensitivity_data <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #34 unique cell lines remain, 7 dups

# save drug sensitivity, cell lines and drugs
gdsc_cl <- colnames(sensitivity_data)
gdsc_dr <- rownames(sensitivity_data)

#################
### CTRP PSet ###
#################

# load in CTRP PSet
load("C:/Users/julia/Desktop/ncRNA Code/data/PSets/CTRP.rds")

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(CTRP, sensitivity.measure = "aac_recomputed",  fill.missing = F))
sensitivity_data <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #32 unique cell lines remain, 7 dups

# save drug sensitivity, cell lines and drugs
ctrp_cl <- colnames(sensitivity_data)
ctrp_dr <- rownames(sensitivity_data)

#################
### CCLE PSet ###
#################

# load in CCLE PSet
ccle <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/CCLE.rds")
ccle <- updateObject(ccle)

# get drug sensitivity data
sensitivity_data <- as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "aac_recomputed",  fill.missing = F))
ccle_sn <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)] #23 unique cell lines remain, 7 dups

# save drug sensitivity, cell lines and drugs
ccle_cl <- colnames(ccle_sn)
ccle_dr <- rownames(ccle_sn)

save(uhnbreast_cl, uhnbreast_dr, uhnbreast2_cl, uhnbreast2_dr, gray_cl, gray_dr, 
        gcsi_cl, gcsi_dr, gdsc_cl, gdsc_dr, ctrp_cl, ctrp_dr, 
        ccle_cl, ccle_dr, file = "temp.RData")

df <- data.frame(PSet = c("UHNBreast", "UHNBreast2", "GRAY", "gCSI", "GDSC2", "CTRP", "CCLE"),
                CL_Count = c(length(uhnbreast_cl), length(uhnbreast2_cl), length(gray_cl), length(gcsi_cl),
                            length(gdsc_cl), length(ctrp_cl), length(ccle_cl)),
                Drug_Count = c(length(uhnbreast_dr), length(uhnbreast2_dr), length(gray_dr), length(gcsi_dr),
                            length(gdsc_dr), length(ctrp_dr), length(ccle_dr)))

#############
### Plots ###
#############

library(reshape2)

### Heat map of cell lines in PSets
uhnbreast2_cl <- unique(samples$sample) #look at code above, samples are all present but mislabeled
all_cl <- unique(c(uhnbreast_cl, uhnbreast2_cl, gray_cl, gcsi_cl, gdsc_cl, ctrp_cl, ccle_cl))

toPlot <- data.frame(
    ubr1 = ifelse(all_cl %in% uhnbreast_cl, 1, 0),
    ubr2 = ifelse(all_cl %in% uhnbreast2_cl, 1, 0),
    gray = ifelse(all_cl %in% gray_cl, 1, 0),
    gcsi = ifelse(all_cl %in% gcsi_cl, 1, 0),
    gdsc = ifelse(all_cl %in% gdsc_cl, 1, 0),
    ctrp = ifelse(all_cl %in% ctrp_cl, 1, 0),
    ccle = ifelse(all_cl %in% ccle_cl, 1, 0))
toPlot$sample <- all_cl
toPlot <- melt(toPlot)

PSet_labels <- c("UHNBreast", "UHNBreast2", "GRAY", "gCSI", "GDSC2", "CTRP", "CCLE")

png("DrugResponse/results/PSetoverlap.png", width = 4, height = 8, res = 600, units = "in")
ggplot(toPlot, aes(x = variable, y = reorder(sample, -as.numeric(factor(sample))), fill = factor(value))) + 
    geom_tile(color = "black") + scale_fill_manual(values = c("0" = "white", "1" = wes_palette("Cavalcanti1")[4]), labels = c("False", "True")) +
    theme_minimal() + scale_x_discrete(labels = PSet_labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "PSet", y = "Breast Cancer Cell Line", fill = "In PSet")
dev.off()



### Heat map of drugs in PSets
all_dr <- unique(c(uhnbreast_dr, uhnbreast2_dr, gray_dr, gcsi_dr, gdsc_dr, ctrp_dr, ccle_dr))

toPlot <- data.frame(
    ubr1 = ifelse(all_dr %in% uhnbreast_dr, 1, 0),
    ubr2 = ifelse(all_dr %in% uhnbreast2_dr, 1, 0),
    gray = ifelse(all_dr %in% gray_dr, 1, 0),
    gcsi = ifelse(all_dr %in% gcsi_dr, 1, 0),
    gdsc = ifelse(all_dr %in% gdsc_dr, 1, 0),
    ctrp = ifelse(all_dr %in% ctrp_dr, 1, 0),
    ccle = ifelse(all_dr %in% ccle_dr, 1, 0))
toPlot$drug <- all_dr
toPlot <- melt(toPlot)

PSet_labels <- c("UHNBreast", "UHNBreast2", "GRAY", "gCSI", "GDSC2", "CTRP", "CCLE")

png("DrugResponse/results/PSetoverlap_drug.png", width = 4, height = 8, res = 600, units = "in")
ggplot(toPlot, aes(x = variable, y = reorder(drug, -as.numeric(factor(drug))), fill = factor(value))) + 
    geom_tile(color = "black") + scale_fill_manual(values = c("0" = "white", "1" = wes_palette("Cavalcanti1")[4]), labels = c("False", "True")) +
    theme_void()
    #theme_minimal() + scale_x_discrete(labels = PSet_labels) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "PSet", y = "Breast Cancer Cell Line", fill = "In PSet")
dev.off()