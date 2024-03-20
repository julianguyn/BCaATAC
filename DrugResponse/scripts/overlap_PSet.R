### This script checks or overlapping cell lines and drugs across PSets
### Also correlations between PSet drug response for overlapping cell lines and drugs

setwd("C:/Users/julia/Documents/BCaATAC")

library(PharmacoGx)
library(ggplot2)
library(wesanderson)
library(reshape2)

# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")


### UHN Breast PSet ###

# load in UHN Breast PSet
uhnbreast <- readRDS("DrugResponse/data/PSet_UHNBreast.rds")
uhnbreast <- updateObject(uhnbreast)

# get drug sensitivity data
ubr1_sen <- as.data.frame(summarizeSensitivityProfiles(uhnbreast, sensitivity.measure = "aac_recomputed",  fill.missing = F))
rm(uhnbreast)


### UHN Breast PSet2 ###

# load in UHN Breast PSet
uhnbreast2 <- readRDS("DrugResponse/data/PharmacoSet.RDS")

# get drug sensitivity data
ubr2_sen <- dcast(uhnbreast2@treatmentResponse$profiles, treatmentid ~ sampleid, value.var = "aac_recomputed")
rownames(ubr2_sen) <- ubr2_sen$treatmentid
ubr2_sen <- ubr2_sen[,-c(1)]
# mapping cell line names
mapping <- c("BT20" = "BT-20", "BT474" = "BT-474", "BT549" = "BT-549", "CAL120" = "CAL-120", "CAL148" = "CAL-148", 
             "CAL51" = "CAL-51", "CAMA1" = "CAMA-1", "EFM19" = "EFM-19", "HDQP1" = "HDQ-P1", "HS578T" = "Hs 578T",
             "JIMT1" = "JIMT-1", "KPL1" = "KPL-1", "LY2" = "MCF-7/LY2", "MCF7" = "MCF-7", "MDAMB157" = "MDA-MB-157",
             "MDAMB175VII" = "MDA-MB-175-VII", "MDAMB231" = "MDA-MB-231", "MDAMB361" = "MDA-MB-361", 
             "MDAMB436" = "MDA-MB-436", "MDAMB468" = "MDA-MB-468", "MX1" = "MX-1", "SKBR3" = "SK-BR-3", 
             "SKBR5" = "SK-BR-5", "SKBR7" = "SK-BR-7", "SUM149" = "SUM149PT", "SUM159" = "SUM159PT", "T47D" = "T-47D",
             "UACC812" = "UACC-812", "ZR751" = "ZR-75-1")
for (i in 1:length(colnames(ubr2_sen))) {
    cell = colnames(ubr2_sen)[i]
    if (cell %in% names(mapping)) {colnames(ubr2_sen)[i] <- unname(mapping[cell])}
}
rm(uhnbreast2)


### Gray PSet ###

# load in Gray PSet
gray <- readRDS("DrugResponse/data/PSet_GRAY2017.rds")
gray <- updateObject(gray)

# get drug sensitivity data
gray_sen <- as.data.frame(summarizeSensitivityProfiles(gray, sensitivity.measure = "aac_recomputed",  fill.missing = F))
rm(gray)


### gCSI PSet ###

# load in gCSI PSet
gcsi <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/gCSI.rds")
gcsi <- updateObject(gcsi)

# get drug sensitivity data
gcsi_sen <- as.data.frame(summarizeSensitivityProfiles(gcsi, sensitivity.measure = "aac_recomputed",  fill.missing = F))
rm(gcsi)

### GDSC PSet ###

# load in GDSC PSet
gdsc <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/GDSC2-8.2.rds")
gdsc <- updateObject(gdsc)

# get drug sensitivity data
gdsc_sen <- as.data.frame(summarizeSensitivityProfiles(gdsc, sensitivity.measure = "aac_recomputed",  fill.missing = F))
rm(gdsc)

### CTRP PSet ###

# load in CTRP PSet
load("C:/Users/julia/Desktop/ncRNA Code/data/PSets/CTRP.rds")

# get drug sensitivity data
ctrp_sen <- as.data.frame(summarizeSensitivityProfiles(CTRP, sensitivity.measure = "aac_recomputed",  fill.missing = F))
rm(CTRP)

### CCLE PSet ###

# load in CCLE PSet
ccle <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/CCLE.rds")
ccle <- updateObject(ccle)

# get drug sensitivity data
ccle_sen <- as.data.frame(summarizeSensitivityProfiles(ccle, sensitivity.measure = "aac_recomputed",  fill.missing = F))
rm(ccle)

###### save sensitivity data ####

# UHN Breast2 Drug Mapping
mapping <- c("carboplatin" = "Carboplatinum", "dipyridamole" = "Dipyridamole", 
             "epirubicin" = "Epirubicin", "eribulin" = "Eribulin", "fluvastatin" = "Fluvastatin", 
             "herceptin" = "Herceptin", "lapatinib" = "Lapatinib", "paclitaxel" = "Paclitaxel")
for (i in 1:length(rownames(ubr2_sen))) {
    cell = rownames(ubr2_sen)[i]
    if (cell %in% names(mapping)) {rownames(ubr2_sen)[i] <- unname(mapping[cell])}
}
save(ubr1_sen, ubr2_sen, gray_sen, gcsi_sen, gdsc_sen, ctrp_sen, ccle_sen, file = "DrugResponse/results/data/sensitivity_data.RData")




#############
### Plots ###
#############

### Drug response correlations between PSets

# function to format drug sensitivity data for correlation
format_sen <- function(pset) {

    # keep only cell lines in analysis
    pset <- pset[,which(colnames(pset) %in% samples$sample)]
    
    # melt dataframe
    pset$drug <- rownames(pset)
    pset <- melt(pset)
    pset$pairs <- paste0(pset$variable, "_", pset$drug)

    return(pset)
}

# function to create a correlation plot and compute correlation coefficient between two psets
corr_pset <- function(pset1, pset2, lab1, lab2) {
    ################################################
    # Inputs:
    #   pset1: drug sensitivity data of first pset
    #   pset2: drug sensitivity data of second pset
    #   lab1: label of first pset
    #   lab2: label of second pset
    # Outputs:
    #   p: correlation graph 
    #   corr_df: will update this dataframe with 
    #            Pearson's correlation coefficient
    #            and number of overlapping features
    ################################################

    # intersected drug-cell pairs
    pset1 <- format_sen(pset1)
    pset2 <- format_sen(pset2)
    
    overlap <- intersect(pset1$pairs, pset2$pairs)
    
    # save intersected drug-cell pairs and order them
    pset1_intersect <- pset1[pset1$pairs %in% overlap,]
    pset1_intersect <- pset1[match(overlap, pset1$pairs),]
    pset2_intersect <- pset2[pset2$pairs %in% overlap,]
    pset2_intersect <- pset2[match(overlap, pset2$pairs),]

    # scatter plot of drug response difference for CTRP and GDSC 
    toPlot <- data.frame(pair = overlap, pset1 = pset1_intersect$value, pset2 = pset2_intersect$value)
    
    # pearson correlation coefficient
    corr <- cor(toPlot$pset1, toPlot$pset2,  method = "pearson", use = "complete.obs")

    p <- ggplot(toPlot, aes(x = pset1, y = pset2)) + 
            geom_smooth(method=lm, show.legend = FALSE, color = "#046C9A") + geom_point(shape = 21, size = 2.5, color = "black", fill = "#899DA4") + 
            geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") + 
            theme_classic() + xlim(c(0, 1)) + ylim(c(0, 1)) + theme(legend.key.size = unit(0.5, 'cm')) +
            labs(x = lab1, y = lab2) + geom_text(x = 0.01, y = 1, label = paste("Corr: ", round(corr, digits = 3)), color = "black")

    print(paste0("Correlation between ", lab1, " and ", lab2, ": ", corr))

    return(p)
}

## TODO: Make a dataframe that stores for each pset pair, number of overlap and correlation

# UHN Breast1
p1 <- corr_pset(ubr1_sen, ubr2_sen/100, "UHNBreast1", "UHNBreast2") #0.813656651551272
p2 <- corr_pset(ubr1_sen, gray_sen, "UHNBreast1", "GRAY") #0.28751603589755
p3 <- corr_pset(ubr1_sen, gcsi_sen, "UHNBreast1", "gCSI") #0.376755118077488
p4 <- corr_pset(ubr1_sen, gdsc_sen, "UHNBreast1", "GDSC2") #0.454090310308515
p5 <- corr_pset(ubr1_sen, ctrp_sen, "UHNBreast1", "CTRP") #0.33775889496196
p6 <- corr_pset(ubr1_sen, ccle_sen, "UHNBreast1", "CCLE") #0.132042976073781

# UHN Breast2
p7 <- corr_pset(ubr2_sen/100, gray_sen, "UHNBreast2", "GRAY") #0.33933595599161
p8 <- corr_pset(ubr2_sen/100, gcsi_sen, "UHNBreast2", "gCSI") #0.516625128943435
p9 <- corr_pset(ubr2_sen/100, gdsc_sen, "UHNBreast2", "GDSC2") #0.346082745902764
p10 <- corr_pset(ubr2_sen/100, ctrp_sen, "UHNBreast2", "CTRP") #0.414636598545278
p11 <- corr_pset(ubr2_sen/100, ccle_sen, "UHNBreast2", "CCLE") #0.473716534195254

# GRAY
p12 <- corr_pset(gray_sen, gcsi_sen, "GRAY", "gCSI") #0.71224637722019
p13 <- corr_pset(gray_sen, gdsc_sen, "GRAY", "GDSC2") #0.450278493429463
p14 <- corr_pset(gray_sen, ctrp_sen, "GRAY", "CTRP") #0.555752295602104
p15 <- corr_pset(gray_sen, ccle_sen, "GRAY", "CCLE") #0.85899408148493

# gCSI
p16 <- corr_pset(gcsi_sen, gdsc_sen, "gCSI", "GDSC2") #0.398261023284974
p17 <- corr_pset(gcsi_sen, ctrp_sen, "gCSI", "CTRP") #0.669311974705811
p18 <- corr_pset(gcsi_sen, ccle_sen, "gCSI", "CCLE") #0.856997604962811

# GDSC
p19 <- corr_pset(gdsc_sen, ctrp_sen, "GDSC2", "CTRP") #0.595141897419167
p20 <- corr_pset(gdsc_sen, ccle_sen, "GDSC2", "CCLE") #0.131232053030648

# CTRP
p21 <- corr_pset(ctrp_sen, ccle_sen, "CTRP", "CCLE") #0.808495639267808

suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

png("DrugResponse/results/figures/overlap_PSet/correlations.png", width = 15, height = 13, res = 600, units = "in")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21,
             ncol = 6, nrow = 6,
             layout_matrix = rbind(c(1, NA, NA, NA, NA, NA),
                                   c(2,  7, NA, NA, NA, NA),
                                   c(3,  8, 12, NA, NA, NA),
                                   c(4,  9, 13, 16, NA, NA),
                                   c(5, 10, 14, 17, 19, NA),
                                   c(6, 11, 15, 18, 20, 21)))
dev.off()




### Heat map of cell lines in PSets

df <- data.frame(PSet = c("UHNBreast", "UHNBreast2", "GRAY", "gCSI", "GDSC2", "CTRP", "CCLE"),
                CL_Count = c(length(uhnbreast_cl), length(uhnbreast2_cl), length(gray_cl), length(gcsi_cl),
                            length(gdsc_cl), length(ctrp_cl), length(ccle_cl)),
                Drug_Count = c(length(uhnbreast_dr), length(uhnbreast2_dr), length(gray_dr), length(gcsi_dr),
                            length(gdsc_dr), length(ctrp_dr), length(ccle_dr)))



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