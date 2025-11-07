# load libraries
suppressPackageStartupMessages({
    library(reshape2)
    library(data.table)
    library(ggplot2)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(GenomicRanges)
    library(dplyr)
    library(limma)
    library(DMRcate)
})

source("utils/palettes.R")
source("utils/get_data.R")
source("utils/methylation_analysis.R")
source("utils/plots/methylation_analysis.R")

###########################################################
# Load in data
###########################################################

# get beta values
betas <- readRDS("data/procdata/TCGA/TCGA_betas.RDS")
rownames(betas) <- betas$CpG
betas$CpG <- NULL

# read in meta data file
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
meta$Sample.Name <- gsub("\\.", "-", meta$Sample.Name)
meta <- meta[match(colnames(betas), meta$Sample.Name), ]

# read in tcga clinical data
pheno <- read.table("data/rawdata/tcga/Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", header = TRUE, row.names = 1)

# get 450k probe annotations
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# get cpgs in each ARCHE
arche_cpgs <- get_arche_cpgs(ann450k)

###########################################################
# Remove ctrl and rs probes
###########################################################

betas <- betas[rownames(betas) %in% rownames(ann450k), ]

###########################################################
# Distribution of beta values per ARCHE
###########################################################

# create toPlot
toPlot <- data.frame(matrix(nrow=0, ncol=3))
for (arche in paste0("ARCHE", 1:6)) {
    samples <- meta$Sample.Name[meta$ARCHE == arche]

    # arche-specific cpg sites
    cpgs <- arche_cpgs$CpGs[arche_cpgs$ARCHE == arche]
    arche_bvals <- betas[rownames(betas) %in% cpgs, colnames(betas) %in% samples]

    # all cpgs for arche tumour
    all_bvals <- betas[, colnames(betas) %in% samples]
    
    bvals <- reshape2::melt(rbind(arche_bvals, all_bvals))
    bvals$ARCHE <- arche
    bvals$label <- c(
        rep("ARCHE\nspecific", ncol(arche_bvals)*nrow(arche_bvals)), 
        rep("All CpG\nsites", ncol(all_bvals)*nrow(all_bvals)))
    toPlot <- rbind(toPlot, bvals)
}

# categorize methylation status
toPlot <- toPlot %>%
    mutate(
        status = case_when(
            value <  0.2 ~ "Hypomethylated",
            value > 0.8 ~ "Hypermethylated",
            between(value, 0.2, 0.8) ~ "Intermediate",
            TRUE ~ NA_character_
        )
    )

# get proportions by methylation status
toPlot <- toPlot %>%
    filter(!is.na(status)) %>%
    group_by(ARCHE, label, status) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(ARCHE, label) %>%
    mutate(prop = n / sum(n))
toPlot$status <- factor(toPlot$status, levels = c("Hypomethylated", "Intermediate", "Hypermethylated"))

# plot beta value proportions
plot_bval_arches(toPlot)

###########################################################
# Format data for DMP
###########################################################

# convert beta values to M values
# add 1e-6 to avoid log 0
mvals <- log2((betas + 1e-6) / (1 - betas + 1e-6))

# format pheno data
pheno <- as.data.frame(t(pheno))
rownames(pheno) <- gsub("\\.", "-", rownames(pheno))
pheno <- pheno[match(meta$Sample.Name, rownames(pheno)),]

# get covariates
meta$age <- pheno$years_to_birth
# meta$gender <- pheno$gender        # all female
# meta$stage <- pheno$pathologic_stage  # drop down to 4 and 7
meta$ARCHE <- gsub("ARCHE", "", meta$ARCHE)
meta$ARCHE <- factor(meta$ARCHE, levels = c(1:6))

###########################################################
# Get differentially methylated probes
###########################################################

# fit lm
design <- model.matrix(~ 0 + ARCHE + age, data = meta)
colnames(design) <- make.names(colnames(design))
fit <- lmFit(mvals, design)

# fit contrasts
contrasts <- makeContrasts(
    ARCHE1 = ARCHE1 - (ARCHE2 + ARCHE3 + ARCHE4 + ARCHE5 + ARCHE6)/5,
    ARCHE2 = ARCHE2 - (ARCHE1 + ARCHE3 + ARCHE4 + ARCHE5 + ARCHE6)/5,
    ARCHE3 = ARCHE3 - (ARCHE1 + ARCHE2 + ARCHE4 + ARCHE5 + ARCHE6)/5,
    ARCHE4 = ARCHE4 - (ARCHE1 + ARCHE2 + ARCHE3 + ARCHE5 + ARCHE6)/5,
    ARCHE5 = ARCHE5 - (ARCHE1 + ARCHE2 + ARCHE3 + ARCHE4 + ARCHE6)/5,
    ARCHE6 = ARCHE6 - (ARCHE1 + ARCHE2 + ARCHE3 + ARCHE4 + ARCHE5)/5,
    levels = design
)
fit2 <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit2)

###########################################################
# Save DMP results
###########################################################

# helper function to save results
save_dmp <- function(arche) {

    res <- topTable(fit, coef = arche, number = Inf)

    # save results
    filename <- paste0("data/results/data/2-MolecularSigAnalysis/DMP/", arche, ".csv")
    write.csv(res, file = filename, quote = FALSE, row.names = FALSE)

    # return sig res
    sig_res <- res[which(abs(res$logFC) > 2 & res$adj.P.Val < 0.05), ]
    message(paste(arche, nrow(sig_res)))
    if (nrow(sig_res) == 0) {
        return(NULL)
    } else {
        sig_res$anno <- ann450k$Relation_to_Island[match(rownames(sig_res), rownames(ann450k))]
        sig_res$ARCHE <- arche
        sig_res$label <- ifelse(sig_res$logFC > 0, "Hypermethylated", "Hypomethylated")
        return(sig_res)
    }
}

dmp1 <- save_dmp("ARCHE1")
dmp2 <- save_dmp("ARCHE2") #43
dmp3 <- save_dmp("ARCHE3") #146
dmp4 <- save_dmp("ARCHE4")
dmp5 <- save_dmp("ARCHE5")
dmp6 <- save_dmp("ARCHE6")

###########################################################
# Plot DMP results
###########################################################

# compile dmps
dmps <- rbind(dmp2, dmp3)

# format counts
toPlot <- dmps %>%
    filter(!is.na(label)) %>%
    count(ARCHE, anno, label)
toPlot$ARCHE <- factor(toPlot$ARCHE, levels = paste0("ARCHE", 3:2))
toPlot$anno <- factor(toPlot$anno, levels = names(dmp_pal))

# sum of count of dmps
counts <- dmps %>%
  count(ARCHE, label)

# plot dmps annotated by cpg island label
plot_dmps_arches(toPlot, counts)

###########################################################
# Get differentially methylated probes
###########################################################

###########################################################
# Differentially methylated regions
###########################################################