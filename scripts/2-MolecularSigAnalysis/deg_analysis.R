# load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})

source("utils/plots/deg_analysis.R")
source("utils/get_data.R")
source("utils/palettes.R")

###########################################################
# Load in data
###########################################################

# load in counts matrix
counts <- get_tcga_rna()

# load in metadata file
meta <- read.csv("data/rawdata/tcga/TCGA_sourcefiles.csv")
meta <- meta[match(colnames(counts), meta$Sample.Name),]

###########################################################
# Un-log2 counts
###########################################################

# counts matrix was log2 normalized, undo that
counts <- round(2^counts - 1)

###########################################################
# Function for DEG
###########################################################

run_DEG <- function(counts, meta, arche = NA, subtype = NA) {

  # create colData
  colData <- data.frame(
    SampleID = meta$Sample.Name,
    ARCHE = ifelse(meta$ARCHE == arche, arche, "Other"),
    Subtype = ifelse(meta$Subtype == subtype, subtype, "Other")
  )
  
  if (!is.na(arche)) {  # if running ARCHEx vs all
    colData$ARCHE <- factor(colData$ARCHE, levels = c("Other", arche))
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design= ~ ARCHE) |>
                                suppressMessages()
    label <- paste0("ARCHE_", arche, "_vs_Other")
  } else {            # if running Subtypex vs all
    colData$Subtype <- factor(colData$Subtype, levels = c("Other", subtype))
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design= ~ Subtype) |>
                                suppressMessages()
    label <- paste0("Subtype_", subtype, "_vs_Other")
  }

  # DEG
  dds <- DESeq(dds) |> suppressMessages()
  res <- results(dds, name=label)
  write.csv(res, file = paste0("data/results/data/2-MolecularSigAnalysis/DEG/", label, ".csv"), quote = F, row.names = T)
  return(res)
}

###########################################################
# Run DEG of each ARCHE vs all
###########################################################

a1 <- run_DEG(counts, meta, arche = "ARCHE1")
a2 <- run_DEG(counts, meta, arche = "ARCHE2")
a3 <- run_DEG(counts, meta, arche = "ARCHE3")
a4 <- run_DEG(counts, meta, arche = "ARCHE4")
a5 <- run_DEG(counts, meta, arche = "ARCHE5")
a6 <- run_DEG(counts, meta, arche = "ARCHE6")

###########################################################
# Run DEG of each Subtype vs all
###########################################################

sb <- run_DEG(counts, meta, subtype = "Basal")
sh <- run_DEG(counts, meta, subtype = "Her2")
sla <- run_DEG(counts, meta, subtype = "LumA")
slb <- run_DEG(counts, meta, subtype = "LumB")
sn <- run_DEG(counts, meta, subtype = "Normal")

###########################################################
# Volcano plots
###########################################################

plot_volcano(a1, "ARCHE1")
plot_volcano(a2, "ARCHE2")
plot_volcano(a3, "ARCHE3")
plot_volcano(a4, "ARCHE4")
plot_volcano(a5, "ARCHE5")
plot_volcano(a6, "ARCHE6")

plot_volcano(sb, "Basal")
plot_volcano(sh, "Her2")
plot_volcano(sla, "LumA")
plot_volcano(slb, "LumB")
#plot_volcano(sn, "Normal") # only 4 sig genes, too small sample size

###########################################################
# Plot MYC volcano plots
###########################################################

plot_volcano_MYC(a2, "ARCHE2")
plot_volcano_MYC(a5, "ARCHE5")
plot_volcano_MYC(sb, "Basal")

###########################################################
# MYC expression analysis across all tumours
###########################################################

myc <- as.data.frame(counts[rownames(counts) == "MYC",]) 
myc$Signature <- meta$Signature[match(rownames(myc), meta$Sample.Name)]
myc$Subtype <- meta$Subtype[match(rownames(myc), meta$Sample.Name)]
myc$Sample.Name <- rownames(myc)
colnames(myc)[1] <- "MYC"

# plot MYC expression
plot_MYCexp(myc, "Signature", "all")
plot_MYCexp(myc, "Subtype", "all")

###########################################################
# ARCHE 2 vs ARCHE 5 MYC expression
###########################################################

# myc expression 
myc <- myc[which(myc$Subtype == "Basal" & myc$Signature != "ARCHE3"),]
plot_MYCexp(myc, "Signature", "basal")

# DEG
myc_counts <- counts[,match(rownames(myc), colnames(counts))]
a25 <- run_DEG(myc_counts, myc, arche = "ARCHE2")
plot_volcano(a25, "Basal ARCHE2")
plot_volcano_MYC(a25, "Basal ARCHE2")