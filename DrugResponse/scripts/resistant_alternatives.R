# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
})

###########################################################
# Load in data
###########################################################

# load class A biomarkers
biomarkers <- read.csv("DrugResponse/results/data/ClassA_Biomarkers.csv")
biomarkers <- biomarkers[which((biomarkers$signature == "Signature5" | biomarkers$signature == "Signature4") & biomarkers$type == "Resistance"),]

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
# Define functions for plotting
###########################################################

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

    print(paste("Drug:", drug))

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
    for (j in 1:nrow(all_drug_sen)) { 
        score = signature_scores[rownames(signature_scores) == signature,colnames(signature_scores) == all_drug_sen$Sample[j]]
        if (is.numeric(score)) {
            all_drug_sen$Score[j] <-  score
        } else {
            all_drug_sen$Score[j] <- NA
        }
    }

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
    png(paste0("DrugResponse/results/figures/pacR/", signature, " - ", drug, ".png"), width=160, height=125, units='mm', res = 600, pointsize=80)
    print({ggplot(all_drug_sen, aes(x = Score, y = AAC, fill = PSet)) + geom_point(size = 3, shape = 21) + 
        geom_smooth(method = "lm", se=F, aes(color = PSet)) + 
        scale_fill_manual(values = pal[names(pal) %in% unique(all_drug_sen$PSet)]) + 
        scale_color_manual(values = pal[names(pal) %in% unique(all_drug_sen$PSet)]) +
        guides(color = 'none', fill = guide_legend(override.aes=list(values = pal[names(pal) %in% unique(all_drug_sen$PSet)], linetype = 0))) +
        theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                                plot.title = element_text(hjust = 0.5, size = 16), legend.key.size = unit(0.7, 'cm')) +
        xlim(-x, x) + ylim(0, 1) + labs(x = paste0(signature, " Similarity Score"), title = paste0(signature, " - ", drug)) 
    })
    dev.off()

    return(corr_res)

}


###########################################################
# Plot candidate pac resistant drugs
###########################################################

# set up palette for plotting
pal <- c("UHNBreast1" = "#588B8B", "UHNBreast2" = "#729B79", "GRAY" = "#BEB2C8", 
         "gCSI" = "#D6CA98", "GDSC2" = "#414073", "CCLE" = "#F37748", "CTRP" = "#A53860")



# create dataframe to store results
res <- data.frame(matrix(nrow=0, ncol=3)) 
colnames(res) <- c("Signature_Drug", "PSet", "Correlation")

# loop through to plot all class B biomarkers
for (i in 1:nrow(biomarkers)) {
    res <- plot_signature_AAC(biomarkers$signature[i], biomarkers$drug[i], res)
}