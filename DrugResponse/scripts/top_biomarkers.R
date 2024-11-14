setwd("C:/Users/julia/Documents/BCaATAC")

library(ggplot2)
library(reshape2)

# set up palette for plotting
pal <- c("UHNBreast1" = "#588B8B", "UHNBreast2" = "#729B79", "GRAY" = "#BEB2C8", 
         "gCSI" = "#D6CA98", "GDSC2" = "#414073", "CCLE" = "#F37748", "CTRP" = "#A53860")


# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")

# read in signature scores
signature_scores <- read.table("Signatures/data/bca_sign.Zscore.txt")
# remove duplicates
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

signature_scores <- signature_scores[,which(colnames(signature_scores) %in% samples$file)]
colnames(signature_scores) <- samples[match(colnames(signature_scores), samples$file), ]$sample
rownames(signature_scores) <- paste0("Signature", 1:6)

# load in drug sensitivity data
load("DrugResponse/results/data/sensitivity_data.RData")
ubr1_sen <- ubr1_sen[,which(colnames(ubr1_sen) %in% samples$sample)] 
ubr2_sen <- ubr2_sen[,which(colnames(ubr2_sen) %in% samples$sample)]
gray_sen <- gray_sen[,which(colnames(gray_sen) %in% samples$sample)] #check this one, 41 but indiv_pset said 42
gcsi_sen <- gcsi_sen[,which(colnames(gcsi_sen) %in% samples$sample)] 
gdsc_sen <- gdsc_sen[,which(colnames(gdsc_sen) %in% samples$sample)] 
ccle_sen <- ccle_sen[,which(colnames(ccle_sen) %in% samples$sample)] 
ctrp_sen <- ctrp_sen[,which(colnames(ctrp_sen) %in% samples$sample)] 


# helper function to get drug sensitivity from each pset
get_drug_sen <- function(pset, drug, df, label) {
    # inputs:
    #    pset: sensitivity data of given pset
    #    drug: string of drug to search
    #    df: compiled dataframe to merge results 
    #    label: pset name

    if (drug %in% rownames(pset)) {
        # save drug response
        res <- melt(pset[rownames(pset) == drug,])
        colnames(res) <- c("Sample", "AAC")
        res$PSet <- label
        # merge results with compiled dataframe
        df <- rbind(df, res)
        return(df)
    } else {return(df)}
}


# function to plot association between signature score and AAC for each biomarker
plot_signature_AAC <- function(signature, drug, corr_res) {

    # set up dataframe to store results
    all_drug_sen <- data.frame(matrix(nrow=0, ncol=3))
    colnames(all_drug_sen) <- c("Sample", "AAC", "PSet")

    # subset to keep only drug of interest
    all_drug_sen <- get_drug_sen(ubr1_sen, drug, all_drug_sen, "UHNBreast1")
    all_drug_sen <- get_drug_sen(ubr2_sen, drug, all_drug_sen, "UHNBreast2")
    all_drug_sen <- get_drug_sen(gray_sen, drug, all_drug_sen, "GRAY")
    all_drug_sen <- get_drug_sen(gcsi_sen, drug, all_drug_sen, "gCSI")
    all_drug_sen <- get_drug_sen(gdsc_sen, drug, all_drug_sen, "GDSC2")
    all_drug_sen <- get_drug_sen(ccle_sen, drug, all_drug_sen, "CCLE")
    all_drug_sen <- get_drug_sen(ctrp_sen, drug, all_drug_sen, "CTRP")

    # get signature scores
    all_drug_sen$Score <- NA
    for (i in 1:nrow(all_drug_sen)) { all_drug_sen$Score[i] <-  signature_scores[rownames(signature_scores) == signature,colnames(signature_scores) == all_drug_sen$Sample[i]]}

    # compute correlation between siganture score and aac for each pset
    for (pset in unique(all_drug_sen$PSet)) {
        tmp <- all_drug_sen[all_drug_sen$PSet == pset,]
        corr <- cor(tmp$AAC, tmp$Score, use="complete.obs", method = "pearson")

        corr_res <- rbind(corr_res, data.frame(Signature_Drug = paste0(signature, " - ", drug), PSet = pset, Correlation = corr))
    }

    # pearson's for all
    print(cor(all_drug_sen$AAC, all_drug_sen$Score, use="complete.obs", method = "pearson"))

    # get x axis limits for plotting
    x <- max(max(all_drug_sen$Score), abs(min(all_drug_sen$Score)))

    # scatter plot
    png(paste0("DrugResponse/results/figures/top_biomarkers/", signature, " - ", drug, ".png"), width=160, height=125, units='mm', res = 600, pointsize=80)
    print({ggplot(all_drug_sen, aes(x = Score, y = AAC, fill = PSet)) + geom_point(size = 3, shape = 21) + 
        geom_smooth(method = "lm", se=F, aes(color = PSet)) + 
        scale_fill_manual(values = pal[names(pal) %in% unique(all_drug_sen$PSet)]) + 
        scale_color_manual(values = pal[names(pal) %in% unique(all_drug_sen$PSet)]) +
        guides(color = 'none', fill = guide_legend(override.aes=list(values = pal[names(pal) %in% unique(all_drug_sen$PSet)], linetype = 0))) +
        theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                                plot.title = element_text(hjust = 0.5, size = 16), legend.key.size = unit(0.7, 'cm')) +
        xlim(-x, x) + ylim(0, 1) + labs(x = paste0(signature, " Similarity Score"), title = paste0(signature, " - ", drug)) })
    dev.off()

    return(corr_res)

}

# read in preclinical biomarkers
classA <- read.csv("DrugResponse/results/data/ClassA_Biomarkers.csv")
classB <- read.csv("DrugResponse/results/data/ClassB_Biomarkers.csv")
classC <- read.csv("DrugResponse/results/data/ClassC_Biomarkers.csv")


### Class A Biomarkers ###

# create dataframe to store results
resA <- data.frame(matrix(nrow=0, ncol=3)) 
colnames(resA) <- c("Signature_Drug", "PSet", "Correlation")

# take top and bottom 5
n=5
plotA <- unique(data.frame(Signature = classA$signature[c(1:5, (nrow(classA)-n):nrow(classA))], Drug = classA$drug[c(1:5, (nrow(classA)-n):nrow(classA))]))

# loop through to plot all class B biomarkers
for (i in 1:nrow(plotA)) {
    resA <- plot_signature_AAC(plotA$Signature[i], plotA$Drug[i], resA)
}

write.csv(resA, file = "DrugResponse/results/data/Top_A_Biomarkers.csv", quote = F, row.names = F)



### Class B Biomarkers ###

# create dataframe to store results
resB <- data.frame(matrix(nrow=0, ncol=3)) 
colnames(resB) <- c("Signature_Drug", "PSet", "Correlation")

plotB <- unique(data.frame(Signature = classB$signature, Drug = classB$drug))

# loop through to plot all class B biomarkers
for (i in 1:nrow(plotB)) {
    resB <- plot_signature_AAC(plotB$Signature[i], plotB$Drug[i], resB)
}

write.csv(resB, file = "DrugResponse/results/data/Top_B_Biomarkers.csv", quote = F, row.names = F)


### Class C Biomarkers ###

# create dataframe to store results
resC <- data.frame(matrix(nrow=0, ncol=3)) 
colnames(resC) <- c("Signature_Drug", "PSet", "Correlation")

plotC <- unique(data.frame(Signature = classC$signature, Drug = classC$drug))

# loop through to plot all class B biomarkers
for (i in 1:nrow(plotC)) {
    resC <- plot_signature_AAC(plotC$Signature[i], plotC$Drug[i], resC)
}

write.csv(resC, file = "DrugResponse/results/data/Top_C_Biomarkers.csv", quote = F, row.names = F)


# Plot signature 3 with paclitaxel
tmp <- data.frame(matrix(nrow=0, ncol=3)) 
plot_signature_AAC("Signature3", "Paclitaxel", tmp)


# TEMP CODE: Plot association with TOP1 inhibitors
drugs <- c("Irinotecan", "Topotecan")

# create dataframe to store results
tmp <- data.frame(matrix(nrow=0, ncol=3)) 
colnames(tmp) <- c("Signature_Drug", "PSet", "Correlation")

# take top and bottom 5
toPlot <- data.frame(Signature = rep(paste0("Signature", 1:6), 2), Drug = c(rep(drugs[1], 6), rep(drugs[2], 6)))

# loop through to plot all class B biomarkers
for (i in 1:nrow(toPlot)) {
    tmp <- plot_signature_AAC(toPlot$Signature[i], toPlot$Drug[i], tmp)
}
