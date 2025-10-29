# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(reshape2)
    library(AnnotationDbi) 
    library(org.Hs.eg.db) 
})

source("utils/get_data.R")
source("utils/mappings.R")
source("utils/process_ccls.R")

###########################################################
# Load and save PSet sensitivity data
###########################################################

# get drug sensitivites
ubr1_sen <- get_drugsen("data/rawdata/psets/PSet_UHNBreast.rds")                             #43 CCLs
ubr2_sen <- get_drugsen("data/rawdata/psets/PharmacoSet.RDS", update = FALSE, map = TRUE)    #42 CCLs
gray_sen <- get_drugsen("data/rawdata/psets/PSet_GRAY2017.rds")                              #42 CCLs
gcsi_sen <- get_drugsen("data/rawdata/psets/gCSI.rds")                                       #25 CCLs
gdsc_sen <- get_drugsen("data/rawdata/psets/GDSC2-8.2.rds")                                  #34 CCLs
ctrp_sen <- get_drugsen("data/rawdata/psets/CTRP.rds", load = TRUE)                          #32 CCLs
ccle_sen <- get_drugsen("data/rawdata/psets/CCLE.rds")                                       #23 CCLs

# save drug sensitivities
save(ubr1_sen, ubr2_sen, gray_sen, gcsi_sen, gdsc_sen, ctrp_sen, ccle_sen,
     file = "data/procdata/CCLs/sensitivity_data.RData")

###########################################################
# Process UBR2 gene counts matrix
###########################################################

# get gene counts matrices from UBR2
ubr2 <- readRDS("data/rawdata/psets/PharmacoSet.RDS")
ubr2 <- ubr2@molecularProfiles@ExperimentList$genes_counts
rna <- ubr2@assays@data$expr
colnames(rna) <- ubr2@colData$sampleid

# map cells, missing: "HBL100" "HCC1008" 
colnames(rna) <- map_cells(colnames(rna))

# keep only cell lines being used
samples <- get_cells()
rna <- rna[,colnames(rna) %in% samples$sample]
df <- data.frame(Genes = rownames(rna), rna)
colnames(df) <- c("Genes", colnames(rna))
# save
write.table(df, file = "data/procdata/CCLs/rna/UBR2_RNA.tsv",  quote = F, sep = "\t", col.names = T, row.names = F)

###########################################################
# Create UBR2 gene counts matrix metadata
###########################################################

# get gene metadata from UBR2
gene_meta <- ubr2@rowRanges

meta <- data.frame(
    Ensembl = gsub("\\..*", "", gene_meta$gene_id),
    Gene.Symbol = gene_meta$gene_name
)

# get entrez ids
meta$EntrezGene.ID = mapIds(
    org.Hs.eg.db,
    keys=meta$Ensembl, 
    column="ENTREZID",
    keytype="ENSEMBL",
    multiVals="first"
)

# save
write.table(meta, file = "data/procdata/CCLs/rna/UBR2_RNA_meta.tsv",  quote = F, sep = "\t", col.names = T, row.names = F)


###########################################################
# Process UBR2 gene counts matrices from other PSets
###########################################################

get_pset_rna("data/rawdata/psets/PSet_GRAY2017.rds", "data/procdata/CCLs/rna/GRAY_RNA.tsv")
get_pset_rna("data/rawdata/psets/gCSI.rds", "data/procdata/CCLs/rna/gCSI_RNA.tsv")
get_pset_rna("data/rawdata/psets/CCLE.rds", "data/procdata/CCLs/rna/CCLE_RNA.tsv")