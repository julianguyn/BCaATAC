setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
    library(readxl)
    library(reshape2)
    library(ggplot2)
})

###########################################################
# Load in TDXd drug response data
###########################################################

# load in TDXd cell line repsonse data
tdxd <- read_excel("DrugResponse/data/TDXd Cell Line Response Data.xlsx", sheet = 1) |>
    as.data.frame()
tdxd <- tdxd[,colnames(tdxd) %in% c("Cell Line", "Average IC50...4")]
#tdxd <- tdxd[,colnames(tdxd) %in% c("Cell Line", "Average IC50...9")]
colnames(tdxd) <- c("Cell.Line", "Avg.IC50")

# remove no response values
tdxd$Avg.IC50[tdxd$Avg.IC50 == "#DIV/0!"] <- NA


###########################################################
# Load in BCa cell line data
###########################################################

# load in BCa cell lines metadata
samples <- read.csv("DrugResponse/data/cl_samples.csv")

# load in signature scores
signature_scores <- read.table("Signatures/data/bca_sign.Zscore.txt")

# remove duplicates
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

# process signature scores
signature_scores <- signature_scores[,which(colnames(signature_scores) %in% samples$file)]
rownames(signature_scores) <- paste0("Signature", 1:6)


###########################################################
# Keep common samples
###########################################################

# match filenames to cell line names
colnames(signature_scores) <- samples$sample[match(colnames(signature_scores), samples$file)]

# manual mapping of mismatches
samples$mapping <- gsub("-", "", samples$sample)
for (i in seq_along(tdxd$Cell.Line)) {
    if (tdxd$Cell.Line[i] %in% samples$mapping) {
        tdxd$Cell.Line[i] <- samples$sample[samples$mapping == tdxd$Cell.Line[i]]
    }
}
tdxd$Cell.Line[tdxd$Cell.Line == "HS578T"] <- "Hs 578T"
tdxd$Cell.Line[tdxd$Cell.Line == "LY2"] <- "MCF-7/LY2"
tdxd$Cell.Line[tdxd$Cell.Line == "MDA-MB-175VII"] <- "MDA-MB-175-VII"
tdxd$Cell.Line[tdxd$Cell.Line == "MPE600"] <- "600MPE"

# keep common samples
common <- intersect(tdxd$Cell.Line, colnames(signature_scores))
tdxd <- tdxd[match(common, tdxd$Cell.Line),]
signature_scores <- signature_scores[,match(common, colnames(signature_scores))]


###########################################################
# Format TDXd response data for CI
###########################################################

rownames(tdxd) <- tdxd$Cell.Line
tdxd$Cell.Line <- NULL
tdxd <- t(tdxd) |> as.data.frame()
rownames(tdxd) <- "TDXd"


###########################################################
# Compute concordance index
###########################################################

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

    # format dataframe for plotting (order by ci and add variable rank)
    combinations <- combinations[order(combinations$ci),]
    combinations$rank <- 1:nrow(combinations)
    combinations$pairs <- paste0(combinations$signature, "-", combinations$drug)
    combinations$pset <- c(rep(label, nrow(combinations)))
    
    return(combinations)
}

tdxd_com <- computeCI(signature_scores, tdxd, "TDXd")


###########################################################
# Plot association for each signature and TDXd
###########################################################

# combine drug response and signature scores
combined <- t(rbind(tdxd, signature_scores)) |> as.data.frame()
combined <- sapply(combined, as.numeric)
combined <- as.data.frame(combined)

# list signatures
signatures <- paste0("Signature", 1:6)

# get axis limits for plotting
x <- max(max(combined[,-c(1)]), abs(min(combined[,-c(1)])))
y <- max(max(combined$TDXd, na.rm=T), abs(min(combined$TDXd, na.rm=T)))

for (signature in signatures) {

    # compute pearson's correlation between signature score and ic50
    corr <- cor(combined$TDXd, combined[,colnames(combined) == signature], 
              use="complete.obs", method = "pearson")

    # create plot for plotting
    toPlot <- combined[,colnames(combined) %in% c("TDXd", signature)]
    colnames(toPlot) <- c("TDXd", "Signature")

    # plot correlation
    png(paste0("DrugResponse/results/figures/tdxd/", signature, ".png"), width=160, height=125, units='mm', res = 600, pointsize=80)
    print({ggplot(toPlot, aes(x = Signature, y = TDXd)) + 
        geom_point(size = 3, shape = 21, fill = "#046C9A") + geom_smooth(method = "lm", se=F, color = "#046C9A") + 
        theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                                plot.title = element_text(hjust = 0.5, size = 16), legend.key.size = unit(0.7, 'cm')) +
        #xlim(-x, x) + ylim(-y, y) + 
        labs(x = paste0(signature, " Similarity Score"), y = "TDXd IC50 Value", 
        title = paste0(signature, "; Pearson's: ", round(corr, 2))) 
    })
    dev.off()
}

