# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
})


###########################################################
# Load in data
###########################################################

# load in drug sensitivity data from PSets
load("DrugResponse/results/data/sensitivity_data.RData")

# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")
samples$file <- gsub("\\.$", "", samples$file)
samples$file <- gsub("\\.", "-", samples$file)

# remove duplicates
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

# load in metadata
meta <- read.csv("MetaData/Lupien/BCa_samples.csv")
meta <- meta[meta$sample %in% samples$file,]


###########################################################
# Get intrinsic molecular subtype
###########################################################

samples$subtype <- meta[match(samples$file, meta$sample),]$subtype

# from https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2911-z/tables/3, set 600MPE to LuminalA
samples[samples$subtype == "cell_line",]$subtype <- "LumA"


###########################################################
# Functions for drug response associations
###########################################################

# function to return dataframe of drug sensitivity for list of drugs and cell lines needed
subset_sen <- function(sen, need_cl, drugs, df, label) {
    sen <- sen[rownames(sen) %in% drugs,]

    if (nrow(sen) > 0) {
        sen$drug <- rownames(sen)
        res <- melt(sen)
        res$label <- label
        res$subtype <- 0
        for (i in 1:nrow(res)) {ifelse(res$variable[i] %in% need_cl, res$subtype[i] <- "Interest", res$subtype[i] <- "NonInterest")}
        df <- rbind(df, res)
    }
    
    return(df)
}

# function calls above on all PSets
get_drugresponse <- function(subtype, drugs) {
    # inputs:
    #   subtype: vector of subtypes from samples$subtype to subset for
    #   drugs: vector of drugs to subset for

    # subset for cell lines
    need_cl <- samples[samples$subtype %in% subtype,]$sample

    # save results in dataframe
    df <- data.frame(matrix(nrow=0, ncol=5))

    # for each pset, subset for the needed cell lines and drugs and save results
    df <- subset_sen(ubr1_sen,need_cl,drugs,df,"UHNBreast1")
    df <- subset_sen(ubr2_sen/100,need_cl,drugs,df,"UHNBreast2")
    df <- subset_sen(gray_sen,need_cl,drugs,df,"GRAY")
    df <- subset_sen(gcsi_sen,need_cl,drugs,df,"gCSI")
    df <- subset_sen(gdsc_sen,need_cl,drugs,df,"GDSC2")
    df <- subset_sen(ctrp_sen,need_cl,drugs,df,"CTRP")
    df <- subset_sen(ccle_sen,need_cl,drugs,df,"CCLE")
    
    #print(table(df$label, df$drug))
    df <- na.omit(df)

    df$label <- factor(df$label, levels = c("UHNBreast1", "UHNBreast2", "GRAY", "gCSI", "GDSC2", "CTRP", "CCLE"))
    return(df)
}


###########################################################
# Plot HER2 drug response assocaitions
###########################################################

drugs <- c("Trastuzumab", "Lapatinib")
df <- get_drugresponse(c("Her2"), drugs)


png("DataExploration/results/figures/Her2-Lapatinib.png", width = 5, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Lapatinib",], aes(x = label, y = value)) +
  geom_boxplot(aes(fill = subtype)) + 
  scale_fill_manual("Subtype", values = c("#ABC8C0", "#9F9F92"), labels = c("Her2", "Other")) +
  theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "PSet", y = "AAC", title = "Her2 and Lapatinib")
dev.off()

png("DataExploration/results/figures/Her2-Trastuzumab.png", width = 3, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Trastuzumab",], aes(x = label, y = value)) +
  geom_boxplot(aes(fill = subtype)) + 
  scale_fill_manual("Subtype", values = c("#ABC8C0", "#9F9F92"), labels = c("Her2", "Other")) +
  theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "PSet", y = "AAC", title = "Her2 and Trastuzumab")
dev.off()


###########################################################
# Plot ER+ drug response associations
###########################################################

drugs <- c("Tamoxifen")
df <- get_drugresponse(c("LumA", "LumB", "Normal"), drugs)

png("DataExploration/results/figures/ER-Tamoxifen.png", width = 3.5, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Tamoxifen",], aes(x = label, y = value)) +
  geom_boxplot(aes(fill = subtype)) + 
  scale_fill_manual("Subtype", values = c("#ABC8C0", "#9F9F92"), labels = c("ER+", "Other")) +
  theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "PSet", y = "AAC", title = "ER+ and Tamoxifen")
dev.off()


###########################################################
# Plot TNBC drug response associations
###########################################################

drugs <- c("Paclitaxel", "Docetaxel", "Doxorubicin", "Epirubicin")
df <- get_drugresponse(c("Basal"), drugs)

png("DataExploration/results/figures/Basal-Docetaxel.png", width = 4, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Docetaxel",], aes(x = label, y = value)) +
  geom_boxplot(aes(fill = subtype)) + 
  scale_fill_manual("Subtype", values = c("#ABC8C0", "#9F9F92"), labels = c("Basal", "Other")) +
  theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "PSet", y = "AAC", title = "Basal and Docetaxel")
dev.off()

png("DataExploration/results/figures/Basal-Doxorubicin.png", width = 3.5, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Doxorubicin",], aes(x = label, y = value)) +
    geom_boxplot(aes(fill = subtype)) + 
    scale_fill_manual("Subtype", values = c("#ABC8C0", "#9F9F92"), labels = c("Basal", "Other")) +
    theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
    labs(x = "PSet", y = "AAC", title = "Basal and Doxorubicin")
dev.off()

png("DataExploration/results/figures/Basal-Epirubicin.png", width = 3.5, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Epirubicin",], aes(x = label, y = value)) +
  geom_boxplot(aes(fill = subtype)) + 
  scale_fill_manual("Subtype", values = c("#ABC8C0", "#9F9F92"), labels = c("Basal", "Other")) +
  theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "PSet", y = "AAC", title = "Basal and Epirubicin")
dev.off()

png("DataExploration/results/figures/Basal-Paclitaxel.png", width = 5, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Paclitaxel",], aes(x = label, y = value)) +
  geom_boxplot(aes(fill = subtype)) + 
  scale_fill_manual("Subtype", values = c("#ABC8C0", "#9F9F92"), labels = c("Basal", "Other")) +
  theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "PSet", y = "AAC", title = "Basal and Paclitaxel")
dev.off()