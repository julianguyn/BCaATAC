# script to assess the predictive power of PAM50


# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(reshape2)
    library(plyr)
    library(genefu)
    library(PharmacoGx)
    library(ggplot2)
    library(ggpubr)
})

############################################################
# Load in data 
############################################################

# load in pam50 scores from DataExploration/pam50_subtypes.R
load("DrugResponse/results/pam50_scores.RData")

# load in drug sensitivity data
load("DrugResponse/results/data/sensitivity_data.RData")


############################################################
# Compute associations with paclitaxel 
############################################################

# set up palette for plotting
pal <- c("UHNBreast1" = "#588B8B", "UHNBreast2" = "#729B79", "GRAY" = "#BEB2C8", 
         "gCSI" = "#D6CA98", "GDSC2" = "#414073", "CCLE" = "#F37748", "CTRP" = "#A53860")


# function to get drug sensitivity from each pset
get_drug_sen <- function(pset, drug, pam50, label) {
    # inputs:
    #    pset: sensitivity data of given pset
    #    drug: string of drug to search
    #    label: pset name
    #    pam50: pam50 scores

    # filter for drug of interest
    res <- melt(pset[rownames(pset) == drug,])
    colnames(res) <- c("Sample", "AAC")

    # merge PAM50 and AAC
    pam50$AAC <- NA
    for (i in 1:nrow(pam50)) { 
        if (rownames(pam50)[i] %in% res$Sample) {
            pam50$AAC[i] <-  res[res$Sample == rownames(pam50)[i],]$AAC
        } else {
            pam50$AAC[i] <- NA
        }
    }
    pam50 <- pam50[complete.cases(pam50$AAC),]
    pam50$PSet <- label
    rownames(pam50) <- NULL
    return(pam50)
}


# function to plot association between PAM50 score and AAC for each biomarker
plot_pam50_AAC <- function(subtype, drug) {

    # subset to keep only drug of interest
    all_drug_sen <- rbind(get_drug_sen(ubr2_sen, drug, gcsi, "UHNBreast2"),
                          get_drug_sen(gray_sen, drug, gray, "GRAY"),
                          get_drug_sen(gcsi_sen, drug, gcsi, "gCSI"),
                          get_drug_sen(ccle_sen, drug, ccle, "CCLE"),
                          get_drug_sen(ctrp_sen, drug, ccle, "CTRP"))

    # keep only subtype of interest
    all_drug_sen <- all_drug_sen[,colnames(all_drug_sen) %in% c(subtype, "AAC", "PSet")]
    colnames(all_drug_sen)[1] <- "Score"

    # create dataframe to store results
    corr_res <- data.frame(matrix(nrow=0, ncol=3)) 
    colnames(corr_res) <- c("Signature_Drug", "PSet", "Correlation")

    # compute correlation between PAM50 subtype and aac for each pset
    for (pset in unique(all_drug_sen$PSet)) {
        tmp <- all_drug_sen[all_drug_sen$PSet == pset,]
        corr <- cor(tmp$AAC, tmp$Score, use="complete.obs", method = "pearson")

        corr_res <- rbind(corr_res, data.frame(Subtype_Drug = paste0(subtype, " - ", drug), PSet = pset, Correlation = corr))
    }

    # pearson's for all
    print(cor(all_drug_sen$AAC, all_drug_sen$Score, use="complete.obs", method = "pearson"))

    # scatter plot
    png(paste0("DrugResponse/results/figures/pam50/", subtype, " - ", drug, ".png"), width=160, height=125, units='mm', res = 600, pointsize=80)
    print({ggplot(all_drug_sen, aes(x = Score, y = AAC, fill = PSet)) + geom_point(size = 3, shape = 21) + 
        geom_smooth(method = "lm", se=F, aes(color = PSet)) + 
        scale_fill_manual(values = pal[names(pal) %in% unique(all_drug_sen$PSet)]) + 
        scale_color_manual(values = pal[names(pal) %in% unique(all_drug_sen$PSet)]) +
        guides(color = 'none', fill = guide_legend(override.aes=list(values = pal[names(pal) %in% unique(all_drug_sen$PSet)], linetype = 0))) +
        theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                                plot.title = element_text(hjust = 0.5, size = 16), legend.key.size = unit(0.7, 'cm')) +
        xlim(0,1) + ylim(0, 1) + labs(x = paste0(subtype, " Score"), title = paste0(subtype, " - ", drug)) })
    dev.off()

    return(corr_res)

}

basal <- plot_pam50_AAC("Basal", "Paclitaxel")
her2 <- plot_pam50_AAC("Her2", "Paclitaxel")
lumA <- plot_pam50_AAC("LumA", "Paclitaxel")
lumB <- plot_pam50_AAC("LumB", "Paclitaxel")
norm <- plot_pam50_AAC("Normal", "Paclitaxel")