# load libraries
suppressPackageStartupMessages({
    library(PharmacoGx)
    library(reshape2)
    library(AnnotationDbi) 
    library(org.Hs.eg.db) 
})

source("source/get_data.R")
source("source/mappings.R")

###########################################################
# Process UBR2 gene counts matrix
###########################################################

# get gene counts matrices from UBR2
ubr2 <- readRDS("DrugResponse/data/PharmacoSet.RDS")
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
write.table(df, file = "Preprocessing/procdata/CCLs/UBR2_RNA.tsv",  quote = F, sep = "\t", col.names = T, row.names = F)

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
write.table(meta, file = "Preprocessing/procdata/CCLs/UBR2_RNA_meta.tsv",  quote = F, sep = "\t", col.names = T, row.names = F)


###########################################################
# Process UBR2 gene counts matrices from other PSets
###########################################################

get_pset_rna("DrugResponse/data/PSet_GRAY2017.rds", "Preprocessing/procdata/CCLs/GRAY_RNA.tsv")
get_pset_rna("DrugResponse/data/gCSI.rds", "Preprocessing/procdata/CCLs/gCSI_RNA.tsv")
get_pset_rna("DrugResponse/data/CCLE.rds", "Preprocessing/procdata/CCLs/CCLE_RNA.tsv")