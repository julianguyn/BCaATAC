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

# read in subtype information
cl_meta <- read.csv("MetaData/cell_line_subtypes.csv")


###########################################################
# Load BCa RNA-Seq data
###########################################################

# load in RNA-Seq counts
uhnbreast2 <- readRDS("DrugResponse/data/PharmacoSet.RDS")
rna <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@assays@data$expr |>
    as.data.frame()
colnames(rna) <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@colData$sampleid

# load in gene metadata
gene_meta <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@rowRanges |> as.data.frame()

# keep only ERBB2
keep <- gene_meta$gene_id[gene_meta$gene_name == "ERBB2"]
erbb2 <- rna[rownames(rna) == keep,]


###########################################################
# Keep common samples
###########################################################

# match filenames to cell line names
colnames(signature_scores) <- samples$sample[match(colnames(signature_scores), samples$file)]

# manual mapping of mismatches in tdxd response
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
# Compute associations for HER2 vs non-HER2 samples
###########################################################

# stratify samples
her2 <- cl_meta$X[cl_meta$true_subtype == "Her2"]
nonher2 <- cl_meta$X[-which(cl_meta$X %in% her2)]
her2 <- her2[-which(is.na(her2))]

# get mapping
her2 <- samples$sample[match(her2, samples$mapping)]
nonher2 <- samples$sample[match(nonher2, samples$mapping)]

# subset signatures
her2 <- signature_scores[,colnames(signature_scores) %in% her2]
nonher2 <- signature_scores[,colnames(signature_scores) %in% nonher2]

# subset drug response
her2_tdxd <- tdxd[,match(colnames(her2), colnames(tdxd))]
nonher2_tdxd <- tdxd[,match(colnames(nonher2), colnames(tdxd))]

# compute concordance index
her2_com <- computeCI(her2, her2_tdxd, "TDXd")
nonher2_com <- computeCI(nonher2, nonher2_tdxd, "TDXd")


###########################################################
# Combine data for plotting
###########################################################

# function to combine drug response and signature scores
combine_data <- function(signature_scores, tdxd) {
    combined <- t(rbind(tdxd, signature_scores)) |> as.data.frame()
    combined <- sapply(combined, as.numeric)
    combined <- as.data.frame(combined)
    return(combined)
}

combined <- combine_data(signature_scores, tdxd)
her2_combined <- combine_data(her2, her2_tdxd)
nonher2_combined <- combine_data(nonher2, nonher2_tdxd)

# list signatures
signatures <- paste0("Signature", 1:6)

###########################################################
# Plot association for each signature and TDXd
###########################################################

# function to plot pearson's correlation
plot_tdxd <- function(combined, path) {
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
        png(paste0(path, signature, ".png"), width=160, height=125, units='mm', res = 600, pointsize=80)
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
}

plot_tdxd(combined, "DrugResponse/results/figures/tdxd/")
plot_tdxd(her2_combined, "DrugResponse/results/figures/tdxd_her2/")
plot_tdxd(nonher2_combined, "DrugResponse/results/figures/tdxd_nonher2/")


###########################################################
# Format RNA-Seq data and keep intersecting cells
###########################################################

#manual mapping of mismatches in rnaseq counts
colnames(erbb2)[colnames(erbb2) == "LY2"] <- "MCF7/LY2"
colnames(erbb2)[colnames(erbb2) == "SUM52"] <- "SUM52PE"
colnames(erbb2)[colnames(erbb2) == "HS578T"] <- "Hs 578T"
colnames(erbb2)[colnames(erbb2) == "SUM149"] <- "SUM149PT"
colnames(erbb2)[colnames(erbb2) == "SUM159"] <- "SUM159PT"

erbb2 <- erbb2[,colnames(erbb2) %in% samples$mapping]
colnames(erbb2) <- samples$sample[match(colnames(erbb2), samples$mapping)]

# keep common samples
common <- intersect(colnames(tdxd), colnames(erbb2))
tdxd <- tdxd[,match(common, colnames(tdxd)),]
signature_scores <- signature_scores[,match(common, colnames(signature_scores))]
erbb2 <- erbb2[,match(common, colnames(erbb2))]


###########################################################
# Format data for linear modelling
###########################################################

# continous erbb2 expression 
erbb2_exp <- as.numeric(erbb2)

# categorized erbb2 expression
erbb2_cat <- erbb2
q <- quantile(as.numeric(erbb2_cat))
erbb2_cat[erbb2_cat <= q[["25%"]]] <- "Q1"
erbb2_cat[erbb2_cat > q[["25%"]] & erbb2_cat <= q[["50%"]]] <- "Q2"
erbb2_cat[erbb2_cat > q[["50%"]] & erbb2_cat <= q[["75%"]]] <- "Q3"
erbb2_cat[erbb2_cat > q[["75%"]] & erbb2_cat <= q[["100%"]]] <- "Q4"
erbb2_cat <- as.character(erbb2_cat) |> as.factor()

# tdxd response
tdxd <- as.numeric(tdxd)

###########################################################
# Linear regression accounting for ERBB2 expression
###########################################################

# function to run linear regression model
run_lm <- function(signature, erbb2) {

    # extract vectors
    scores <- signature_scores[rownames(signature_scores) == signature,] |>
        as.numeric()

    # run model
    res <- lm(response ~ scores + erbb2)
    print(summary(res))
}

# ERBB2 Continuoue Expression
run_lm("Signature1", erbb2_exp)
run_lm("Signature2", erbb2_exp)
run_lm("Signature3", erbb2_exp)
run_lm("Signature4", erbb2_exp)
run_lm("Signature5", erbb2_exp)
run_lm("Signature6", erbb2_exp)

# ERBB2 Categories
run_lm("Signature1", erbb2_cat)
run_lm("Signature2", erbb2_cat)
run_lm("Signature3", erbb2_cat)
run_lm("Signature4", erbb2_cat)
run_lm("Signature5", erbb2_cat)
run_lm("Signature6", erbb2_cat)


###########################################################
# Plot scatter plot by HER2 expression groups
###########################################################

merged <- cbind(t(signature_scores), data.frame(IC50 = response, Group = erbb2_cat))

# function to create scatter plot by signature
plot_scatter <- function(signature) {

    toPlot <- merged[,colnames(merged) %in% c(signature, "IC50", "Group")]
    colnames(toPlot)[1] <- "Score"

    ggplot(toPlot, aes(x = Score, y = IC50, fill = Group)) + geom_point(size = 3, shape = 21) + 
        geom_smooth(method = "lm", se=F, aes(color = Group)) + 
        guides(color = 'none') +
        theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                                legend.key.size = unit(0.7, 'cm')) +
        labs(x = paste(signature, "Similarity Score"), fill = "ERBB2 Group")
}

plot_scatter("Signature1")
plot_scatter("Signature2")
plot_scatter("Signature3")
plot_scatter("Signature4")
plot_scatter("Signature5")
plot_scatter("Signature6")