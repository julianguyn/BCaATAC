# load libraries
suppressPackageStartupMessages({
  library(maftools)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  library(ggh4x)
  library(ggsignif)
})

source("utils/get_data.R")
source("utils/plots/mutation_analysis.R")
source("utils/bca_mutations.R")
source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# read in meta data file
meta <- read.csv("data/rawdata/TCGA/TCGA_sourcefiles.csv")
meta <- meta[!is.na(meta$SNV.File.Name),]
meta$SNV.File.Name <- paste0("data/rawdata/tcga/maf/", meta$SNV.File.Name)

# load in matrix file from NMF
mat <- get_arche_tcga()
mat <- mat[mat$variable %in% meta$ATAC.Seq.File.Name,]

# get mafs
mafs <- merge_mafs(meta$SNV.File.Name)

###########################################################
# MAF summaries all tumours
###########################################################

plot_mafSummary(arche = "all")
plot_mafSummary(arche = "ARCHE1")
plot_mafSummary(arche = "ARCHE2")
plot_mafSummary(arche = "ARCHE3")
plot_mafSummary(arche = "ARCHE4")
plot_mafSummary(arche = "ARCHE5")
plot_mafSummary(arche = "ARCHE6")

###########################################################
# Compute and plot TMB
###########################################################

# get tmb and format plot
tmb <- tmb(mafs)
tmb$Sample_ID <- sub("^((?:[^-]+-){2}[^-]+).*", "\\1", tmb$Tumor_Sample_Barcode)
tmb$ARCHE <- meta$ARCHE[match(tmb$Sample_ID, gsub("\\.", "-", meta$Sample.Name))]

# ANOVA
tmb.aov <- aov(total_perMB_log ~ ARCHE, data = tmb)
summary(tmb.aov)
#            Df Sum Sq Mean Sq F value Pr(>F)  
#ARCHE        5  1.441 0.28815    3.07 0.0153 *
#Residuals   63  5.913 0.09385                 

# Tukey's HSD
tukey <- TukeyHSD(tmb.aov)$ARCHE |> as.data.frame()
sig_tukey <- tukey[tukey$"p adj" < 0.05,]

# plot TMB ~ ARCHE
plot_tmb_boxplot(tmb)
plot_tmb_waterfall(tmb)

###########################################################
# Create mutations counts matrix
###########################################################

# create counts matrix
mut <- mutCountMatrix(
  mafs,
  includeSyn = FALSE,
  countOnly = NULL,
  removeNonMutated = TRUE
)

write.table(mut, 
            file = "data/results/data/2-MolecularSigAnalysis/TCGA_mutation_matrix.tsv", 
            quote = F, sep = "\t", col.names = T, row.names = T)

###########################################################
# Map files to metadata
###########################################################

mapping <- data.frame(snv = colnames(mut))
mapping$Sample.Name <- sub("^([^-]+-[^-]+-[^-]+).*", "\\1", mapping$snv)
mapping$Sample.Name <- gsub("-", "\\.", mapping$Sample.Name)

meta$snv_label <- mapping$snv[match(meta$Sample.Name, mapping$Sample.Name)]
write.csv(meta, file = "metadata/TCGA_mutation_meta.csv", quote = F, row.names = F)

###########################################################
# Formating dataframe for plotting
###########################################################

# filter for BCa-relevant mutations of interest
toPlot <- mut[rownames(mut) %in% bca_mutations,]

# plot presence of BCa-relevant mutations
plot_BCa_mutations(toPlot)
