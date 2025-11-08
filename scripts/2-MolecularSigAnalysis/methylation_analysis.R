# load libraries
suppressPackageStartupMessages({
    library(reshape2)
    library(data.table)
    library(ggplot2)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(GenomicRanges)
    library(dplyr)
    library(limma)
    library(missMethyl)
    library(DMRcate)
    library(qusage)
    library(ggpubr)
    library(grid)
    library(gridExtra)
    library(ggh4x)
})

source("utils/palettes.R")
source("utils/get_data.R")
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

# create gr of anno
anno_gr <- GRanges(
    seqnames = ann450k$chr,
    ranges = IRanges(start = ann450k$pos, end = ann450k$pos),
    probeID = rownames(ann450k)
)

# get cpgs in each ARCHE
arche_cpgs <- get_arche_cpgs(ann450k, anno_gr)

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

# plot beta value proportions
plot_bval_arches(toPlot)

###########################################################
# Format data for DMP
###########################################################

# convert beta values to M values (add 1e-6 to avoid log 0)
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
# Gene ontology of DMPs
###########################################################

#go2 <- gometh(sig.cpg = rownames(dmp2), all.cpg = rownames(mvals), collection = "GO")
#go3 <- gometh(sig.cpg = rownames(dmp3), all.cpg = rownames(mvals), collection = "GO")
# all FDR = 1

###########################################################
# Get CpG annotations (for dmrcate)
###########################################################

# helper function to get cpg annotation
get_annot <- function(arche) {

    myannotation <- cpg.annotate(
        datatype = "array",
        object = as.matrix(mvals),
        what = "M",
        analysis.type = "differential",
        design = design,
        contrasts = TRUE,
        cont.matrix = contrasts,
        coef = arche,
        arraytype = "450K",
        fdr = 0.1
    )
    return(myannotation)
}

# get cpg annotations
arche2_annot <- get_annot("ARCHE2")
arche3_annot <- get_annot("ARCHE3")

###########################################################
# Get differentially methylated regions
###########################################################

dmr2 <- extractRanges(
    dmrcate(arche2_annot, lambda = 1000, C = 2),    #1026 dmps, 162 DMRs
    genome = "hg19"
)
dmr3 <- extractRanges(
    dmrcate(arche3_annot, lambda = 1000, C = 2),    #409 dmps, 52 DMRs
    genome = "hg19"
)
save(dmr2, dmr3, file = "data/procdata/TCGA/DMRs.RData")
# no DMPs (and hence DMRs) for other ARCHEs

###########################################################
# Quick gene region analysis with missMethyl
###########################################################

# gene ontology
go2 <- goregion(dmr2, all.cpg=rownames(mvals), collection="GO", array.type="450K")
go3 <- goregion(dmr3, all.cpg=rownames(mvals), collection="GO", array.type="450K")
go3 <- go3[go3$FDR < 0.05,] # cell adhesion, 3 terms only, none for a2

# kegg pathway (none with FDR<0.05)
kegg2 <- goregion(dmr2, all.cpg=rownames(mvals), collection="KEGG", array.type="450K")
kegg3 <- goregion(dmr3, all.cpg=rownames(mvals), collection="KEGG", array.type="450K")

# hallmarks (ran this on November 7, 2025 @ 8pm est)
# none with FDR<0.05
hallmark <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds"))
hallmarks2 <- gsaregion(dmr2, all.cpg=rownames(mvals), collection=hallmark)
hallmarks3 <- gsaregion(dmr3, all.cpg=rownames(mvals), collection=hallmark)

# for arche2:                                 N DE      P.DE FDR
#HALLMARK_MYC_TARGETS_V1                    200  0 1.0000000   1
#HALLMARK_MYC_TARGETS_V2                     58  0 1.0000000   1

###########################################################
# Plot DMRs
###########################################################

# get colours for each TCGA samples (by ARCHE)
cols <- ARCHE_pal[as.character(paste0("ARCHE", meta$ARCHE))]

# plot mean diff
plot_all_dmrs(dmr2, "ARCHE2", w = 16, h = 5)
plot_all_dmrs(dmr3, "ARCHE3", w = 8, h = 5)

# plot ARCHE2 and FOXA1
plot_dmr_arche(dmr2, arche2_annot, "ARCHE2", 103, "FOXA1", h = 12)

###########################################################
# Find MYC target genes in ARCHE2
###########################################################

# get myc target genes
myc_targs <- read.gmt("data/rawdata/gmt/All_MYC_Target_Signatures.gmt") 
myc_genes <- unique(unlist(myc_targs))

# get genes in ARCHE2 DMRs
arche2_genes <- dmr2$overlapping.genes[!is.na(dmr2$overlapping.genes)]
arche2_genes <- unique(trimws(unlist(strsplit(arche2_genes, ","))))

# get ARCHE2 DMR genes in myc genes
overlap <- arche2_genes[arche2_genes %in% myc_genes]

# get ARCHE3 genes too while we're at it
arche3_genes <- dmr3$overlapping.genes[!is.na(dmr3$overlapping.genes)]
arche3_genes <- unique(trimws(unlist(strsplit(arche3_genes, ","))))

# save all genes for DEG analysis later
save(arche2_genes, arche3_genes, overlap, file = "data/results/data/2-MolecularSigAnalysis/DMR_genes.RData")

###########################################################
# Plot dmrcate plots of MYC target genes in ARCHE2
###########################################################

for (gene in overlap) {
    plot_dmr_arche(dmr2, arche2_annot, "ARCHE2", grep(gene, dmr2$overlapping.genes), gene, h = 15)
}

###########################################################
# Plot gene overlaps of MYC Gene Sets and ARCHE2 DMRs
###########################################################

plot_arche2_myc_dmrs(overlap, myc_targs)