# script to do sanity checks of correlations between known drug and molecular subtypes

setwd("C:/Users/julia/Documents/BCaATAC")

library(ggplot2)
library(reshape2)

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

# get intrinsic molecular subtype
samples$subtype <- meta[match(samples$file, meta$sample),]$subtype

# from https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2911-z/tables/3, set 600MPE to LuminalA
samples[samples$subtype == "cell_line",]$subtype <- "LumA"


# function to return dataframe of drug sensitivity for list of drugs and cell lines needed
subset_sen <- function(sen, need_cl, drugs, df, label) {
    tmp <- sen[rownames(sen) %in% drugs, colnames(sen) %in% need_cl]
    if (nrow(tmp) > 0) {
        tmp$drug <- rownames(tmp)
        res <- melt(tmp)
        res$label <- label
        df <- rbind(df, res)
    }
    return(df)
}

get_drugresponse <- function(subtype, drugs) {
    # inputs:
    #   subtype: vector of subtypes from samples$subtype to subset for
    #   drugs: vector of drugs to subset for

    # subset for cell lines
    need_cl <- samples[samples$subtype %in% subtype,]$sample

    # save results in dataframe
    df <- data.frame(matrix(nrow=0, ncol=4))

    # for each pset, subset for the needed cell lines and drugs and save results
    df <- subset_sen(ubr1_sen,need_cl,drugs,df,"UHNBreast1")
    df <- subset_sen(ubr2_sen/100,need_cl,drugs,df,"UHNBreast2")
    df <- subset_sen(gray_sen,need_cl,drugs,df,"GRAY")
    df <- subset_sen(gcsi_sen,need_cl,drugs,df,"gCSI")
    df <- subset_sen(gdsc_sen,need_cl,drugs,df,"GDSC2")
    df <- subset_sen(ctrp_sen,need_cl,drugs,df,"CTRP")
    df <- subset_sen(ccle_sen,need_cl,drugs,df,"CCLE")
    
    print(table(df$label, df$drug))
    df <- na.omit(df)

    df$label <- factor(df$label, levels = c("UHNBreast1", "UHNBreast2", "GRAY", "gCSI", "GDSC2", "CTRP", "CCLE"))
    return(df)
}


## HER2
drugs <- c("Trastuzumab", "Lapatinib")
df <- get_drugresponse(c("Her2"), drugs)


png("DataExploration/results/figures/Her2-Lapatinib.png", width = 3, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Lapatinib",], aes(x = label, y = value)) +
    geom_hline(yintercept = 0.5, color = "#ABC8C0", linetype = "dashed") +
    geom_boxplot(fill = "#ABC8C0") + geom_jitter(shape = 16, position=position_jitter(0.2)) +
    theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "PSet", y = "AAC", title = "Her2-Lapatinib")
dev.off()

png("DataExploration/results/figures/Her2-Trastuzumab.png", width = 3, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Trastuzumab",], aes(x = label, y = value)) +
    geom_hline(yintercept = 0.5, color = "#ABC8C0", linetype = "dashed") +
    geom_boxplot(fill = "#ABC8C0") + geom_jitter(shape = 16, position=position_jitter(0.2)) +
    theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "PSet", y = "AAC", title = "Her2-Trastuzumab")
dev.off()


#ER+
drugs <- c("Tamoxifen")
df <- get_drugresponse(c("LumA", "LumB", "Normal"), drugs)

png("DataExploration/results/figures/ER-Tamoxifen.png", width = 3, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Tamoxifen",], aes(x = label, y = value)) +
    geom_hline(yintercept = 0.5, color = "#ABC8C0", linetype = "dashed") +
    geom_boxplot(fill = "#ABC8C0") + geom_jitter(shape = 16, position=position_jitter(0.2)) +
    theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "PSet", y = "AAC", title = "ER-Tamoxifen")
dev.off()

#TNBC
drugs <- c("Paclitaxel", "Docetaxel", "Doxorubicin", "Epirubicin")
df <- get_drugresponse(c("Basal"), drugs)

png("DataExploration/results/figures/Basal-Docetaxel.png", width = 3, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Docetaxel",], aes(x = label, y = value)) +
    geom_hline(yintercept = 0.5, color = "#ABC8C0", linetype = "dashed") +
    geom_boxplot(fill = "#ABC8C0") + geom_jitter(shape = 16, position=position_jitter(0.2)) +
    theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "PSet", y = "AAC", title = "Basal-Docetaxel")
dev.off()

png("DataExploration/results/figures/Basal-Doxorubicin.png", width = 3, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Doxorubicin",], aes(x = label, y = value)) +
    geom_hline(yintercept = 0.5, color = "#ABC8C0", linetype = "dashed") +
    geom_boxplot(fill = "#ABC8C0") + geom_jitter(shape = 16, position=position_jitter(0.2)) +
    theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "PSet", y = "AAC", title = "Basal-Doxorubicin")
dev.off()

png("DataExploration/results/figures/Basal-Epirubicin.png", width = 3, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Epirubicin",], aes(x = label, y = value)) +
    geom_hline(yintercept = 0.5, color = "#ABC8C0", linetype = "dashed") +
    geom_boxplot(fill = "#ABC8C0") + geom_jitter(shape = 16, position=position_jitter(0.2)) +
    theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "PSet", y = "AAC", title = "Basal-Epirubicin")
dev.off()

png("DataExploration/results/figures/Basal-Paclitaxel.png", width = 3, height = 5, res = 600, units = "in")
ggplot(df[df$drug == "Paclitaxel",], aes(x = label, y = value)) +
    geom_hline(yintercept = 0.5, color = "#ABC8C0", linetype = "dashed") +
    geom_boxplot(fill = "#ABC8C0") + geom_jitter(shape = 16, position=position_jitter(0.2)) +
    theme_classic() + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "PSet", y = "AAC", title = "Basal-Paclitaxel")
dev.off()