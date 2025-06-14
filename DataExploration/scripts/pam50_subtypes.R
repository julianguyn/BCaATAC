# Script to call PAM50 molecular subtypes from cell line RNA-Seq data

setwd("C:/Users/julia/Documents/BCaATAC")

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

# load in PAM50
data(pam50.robust)

# get original subtype annotation
subtype_meta <- fread("DataExploration/data/bcacells_annotation.tsv")
subtype_meta <- data.frame(Cell = subtype_meta$sample, Subtype = subtype_meta$subtype)

# load in RNA-Seq metadata (to get gene names)
uhnbreast2 <- readRDS("DrugResponse/data/PharmacoSet.RDS")
gene_meta <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@rowRanges |> as.data.frame()
rm(uhnbreast2)

# format colnames of gene_meta
gene_meta <- gene_meta[,c("gene_id", "gene_name")]
colnames(gene_meta) <- c("GeneID", "Gene.symbol")

# load in cell line rnaseq data
load("DataExploration/data/rnaseq_melted.RData")

# average 'UACC812' replicate in GRAY
rep <- gray[gray$Var2 == 'UACC812',]
rep$value[1:49291] <- (rep$value[1:49291] + rep$value[49292:98582]) / 2
gray <- rbind(gray[-which(gray$Var2 == 'UACC812'),], rep[1:49291,])

# load in GDSC pset
#gdsc <- readRDS("C:/Users/julia/Desktop/ncRNA Code/data/PSets/GDSC2-8.2.rds") |> updateObject()
#gdsc <- summarizeMolecularProfiles(gdsc, mDataType='Kallisto_0.46.1.rnaseq.counts') |> assay() 
#gdsc <- melt(gdsc)
#gdsc <- gdsc[gdsc$Var2 %in% subtype_meta$Cell,]
# too many NA values, unable to compute PAM50

############################################################
# Compute PAM50 subtypes
############################################################

# function to compute subtypes using genefu
pam50subtype <- function(melted_df) {

    # convert dataframe into feature~sample matrix
    df <- reshape2::dcast(melted_df, Var2 ~ Var1)
    rownames(df) <- df$Var2
    df <- df[,-c(1)]

    # match and get gene names
    colnames(df) <- gene_meta$Gene.symbol[match(colnames(df), gene_meta$GeneID)]

    # get subtype predictions for PAM50
    SubtypePredictions <- molecular.subtyping(sbt.model = "pam50", data = df,
                                          annot = gene_meta, do.mapping = FALSE)

    # reutrn subtypes and subtyping probability
    res_df <- cbind(Subtype = as.character(SubtypePredictions$subtype), 
                    as.data.frame(SubtypePredictions$subtype.proba))
    return(res_df)
}

ubr2 <- pam50subtype(cells)
ccle <- pam50subtype(ccle)
gcsi <- pam50subtype(gcsi)
gray <- pam50subtype(gray)

save(ubr2, ccle, gcsi, gray, file = "DrugResponse/results/pam50_scores.RData")

############################################################
# Compile results to plot
############################################################

# compile subtype results
toPlot <- rbind.fill(
    as.data.frame(t(ubr2))[rownames(t(ubr2)) == "Subtype",], 
    as.data.frame(t(gray))[rownames(t(gray)) == "Subtype",],
    as.data.frame(t(gcsi))[rownames(t(gcsi)) == "Subtype",],
    as.data.frame(t(ccle))[rownames(t(ccle)) == "Subtype",]
) |> t() |> as.data.frame()
colnames(toPlot) <- c("UBR2", "GRAY", "gCSI", "CCLE")

# save true PAM50 results
true_subtype <- c("ER+", "Her2", "Her2", "Her2", "TNBC", "TNBC", "ER+",
                  "ER+", "Her2", "TNBC", "Her2", "TNBC", "NA", "ER+",
                  "TNBC", "TNBC", "TNBC", "Her2", "ER+", "ER+", "Her2", 
                  "ER+", "TNBC", "TNBC", "TNBC", "TNBC", "TNBC", "TNBC",
                  "ER+", "TNBC", "Her2", "TNBC", "TNBC", "TNBC", "TNBC",
                  "ER+", "ER+", "TNBC", "TNBC", "TNBC", "TNBC", "TNBC")
# MCF10A is a TN non-tumorigenetic breast cancer cell
#https://pmc.ncbi.nlm.nih.gov/articles/PMC5665029/
toPlot$true_subtype <- true_subtype

# save results
write.csv(toPlot, file = "MetaData/cell_line_subtypes.csv", quote = F, row.names = T)

# format dataframe for plotting
toPlot$samples <- rownames(toPlot)
toPlot <- melt(toPlot, id.vars = "samples")


############################################################
# Plot PAM50 subtypes
############################################################

# function for plotting tiles
plot_tiles <- function(data, x, x_lab, y_lab = TRUE, legend = FALSE) {
    y_title <- ifelse(y_lab == TRUE, "Cell Line", "")
    f_title <- ifelse(legend == FALSE, "", "Subtype")

    p <- ggplot(data, aes(x = x, y = samples, fill = value)) + 
    geom_tile(color = "white") + theme_void() +
    geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_manual(values = subtype_pal, na.value = "#D1D7DD") +        
    theme(
            axis.text.x = element_text(vjust = 0.5),
            axis.title.x = element_text(),
            axis.text.y = element_text(vjust = 0.5, hjust = 1),
            axis.title.y = element_text(angle = 90),
            strip.text.x = element_text(),
            legend.key.size = unit(0.7, 'cm')
        ) +
    labs(x = x_lab, y = y_title, fill = f_title)

    if (y_lab == FALSE) {
        p <- p + theme(axis.text.y = element_blank())
    }
    if (legend == FALSE) {
        p <- p + guides(fill="none") + 
            theme(axis.text.x = element_blank())
    }

    return(p)
}

# set colours for plotting
subtype_pal <- c("Basal" = "#AF4C5B", "TNBC" = "#AF4C5B", 
                 "Her2" = "#EED4D3", "ER+" = "#6567AA", 
                 "LumA" = "#B3B4D0", "LumB" = "#363E62", 
                 "Normal" = "#6365AF", "Not Available" = "#eFeBF7")

# plot true PAM50 subtypes
png("DataExploration/results/figures/pam50_true.png", width = 2, height = 6, res = 600, units = "in")
plot_tiles(
    data = toPlot[toPlot$variable == "true_subtype",], 
    x = 1, 
    x_lab = "\nTrue\nSubtype"
)
dev.off()

# plot PAM50 subtypes
png("DataExploration/results/figures/pam50_pset.png", width = 6, height = 6, res = 600, units = "in")
plot_tiles(
    data = toPlot[toPlot$variable != "true_subtype",], 
    x = toPlot[toPlot$variable != "true_subtype",]$variable, 
    x_lab = "\nPSet", 
    y_lab = FALSE, 
    legend = TRUE
)
dev.off()

# combine figures