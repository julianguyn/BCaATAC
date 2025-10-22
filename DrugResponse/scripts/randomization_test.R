# randomizes the drug sensitivity to see if biomarker detection is better than chance alone

suppressMessages(library(survcomp))
library(wesanderson)
library(ggplot2)
library(ggh4x)
library(reshape2)
library(dplyr)


# load in BCa cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")

# load in sensitivity data
load("DrugResponse/results/data/sensitivity_data.RData")

# get signature scores
signature_scores <- read.table("Signatures/data/bca_sign.Zscore.txt")
# remove duplicates
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

signature_scores <- signature_scores[,which(colnames(signature_scores) %in% samples$file)]

save(ubr1_sen, ubr2_sen, gray_sen, gcsi_sen, gdsc_sen, ctrp_sen, ccle_sen, signature_scores, samples,
     file = "temp_sendtoH4H.RData")


### Function to compute CI ### (modified to return just number of classA biomarkers)
computeCI <- function(signature_scores, samples, sensitivity_data) {

    # subset for overlapping samples in the pset
    pset_sig <- signature_scores[,which(colnames(signature_scores) %in% samples[which(samples$sample %in% colnames(sensitivity_data)),]$file)]
    pset_sen <- sensitivity_data[,which(colnames(sensitivity_data) %in% samples$sample)]

    # shuffle drug sensitivity data (shuffle by drug across the cell lines)
    cell_lines <- colnames(pset_sen)  
    pset_sen <- as.data.frame(t(sapply(as.data.frame(t(pset_sen)), sample)))
    colnames(pset_sen) <- cell_lines

    # create data frame to hold results
    combinations <- as.data.frame(matrix(data = NA, nrow = nrow(pset_sig) * nrow(pset_sen), ncol = 7))
    colnames(combinations) <- c("signature", "drug", "ci", "pvalue", "se", "upper", "lower")
    combinations$signature <- rep(rownames(pset_sig), nrow(pset_sen))
    combinations$drug <- rep(rownames(pset_sen), each = nrow(pset_sig))

    # compute concordance index
    for (i in 1:nrow(combinations)){
        #print(paste0(i, " out of ", nrow(combinations), " complete"))

        ci <- survcomp::concordance.index(as.numeric(pset_sen[combinations[,2][i],]), 
                                            # row from pset_sen with the drug for all samples
                                            surv.time = as.numeric(unlist(-pset_sig[combinations[,1][i],])), 
                                            # row of signature i in signature scores matrix with cell lines as columns
                                            surv.event = rep(1,length(pset_sen)), 
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
   
    num_sig_pairs <- nrow(combinations[(combinations$ci > 0.7 | combinations$ci < 0.3) & combinations$FDR < 0.05,])
    if (class(num_sig_pairs) == "NULL") {
        return(0)
    } else {
        return(num_sig_pairs)
    }
}



# save number of significant transcripts
res <- as.data.frame(matrix(data = NA, nrow = 50, ncol = 7))
colnames(res) <- c("UHNBreast1","UHNBreast2", "GRAY", "gCSI", "GDSC", "CTRP", "CCLE")


for (i in 1:1000) {
  
  if (i %% 10 == 0) {cat(paste("On iteration ",i,"\n",sep=""))}

  ubr1 <- computeCI(signature_scores, samples, ubr1_sen)
  ubr2 <- computeCI(signature_scores, samples, ubr2_sen)
  gray <- computeCI(signature_scores, samples, gray_sen)
  gcsi <- computeCI(signature_scores, samples, gcsi_sen)
  gdsc <- computeCI(signature_scores, samples, gdsc_sen)
  ctrp <- computeCI(signature_scores, samples, ctrp_sen)
  ccle <- computeCI(signature_scores, samples, ccle_sen)
   
  row <- data.frame(ubr1, ubr2, gray, gcsi, gdsc, ctrp, ccle)
  res[i,] <- row

}


save(res, file = "randomize_BCaATAC.RData")

######

# load in randomization results from H4H
load("DrugResponse/results/data/randomize_BCaATAC.RData")
toPlot <- melt(res)

# load in true number of class A biomarkers
classA <- read.csv("DrugResponse/results/data/ClassA_Biomarkers.csv")

# count number of biomarkers for each pset
count_df <- classA %>% group_by(pset) %>% summarise(count = n())
count_df$label <- c("CCLE", "CTRP", "GDSC", "GRAY", "UHNBreast1")



# function to create histogram of distribution of randomization results
plot_randomization <- function(pset) {

    # subset df to only pset of interest
    df <- toPlot[toPlot$variable == pset,]

    # get true number of associations
    true <- count_df[count_df$label == pset,]$count

    # compute p-value
    pval <- nrow(df[df$value >= true,]) / nrow(df)

    p <- ggplot(df, aes(x = value)) + 
        geom_histogram(binwidth = 1, color = "black", fill = "#899DA4") +
        geom_vline(xintercept = true, linetype = "dashed") +
        annotate("text", x = Inf, y = Inf, label = paste0("P-Value: ", pval), hjust = 1.1, vjust = 1.5, size = 5) +
        theme_classic() + labs(x = "Count of Class A Associations", y = "Density", title = pset)

    return(p)

}

p1 <- plot_randomization("UHNBreast1")
p2 <- plot_randomization("GRAY")
p3 <- plot_randomization("GDSC")
p4 <- plot_randomization("CTRP")
p5 <- plot_randomization("CCLE")


suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

png("DrugResponse/results/figures/randomization_check.png", width = 12, height = 8, res = 600, units = "in")
grid.arrange(p1, p2, p3, p4, p5,
             ncol = 3, nrow = 2,
             layout_matrix = rbind(c(1, 2, 3),
                                   c(4,  5, NA)))
dev.off()

