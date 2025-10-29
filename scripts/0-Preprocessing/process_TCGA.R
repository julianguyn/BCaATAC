# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
})


###########################################################
# Load in metadata
###########################################################

# load in tcga metadata
meta <- read_excel("data/rawdata/tcga/TCGA-ATAC_DataS1_DonorsAndStats_v4.xlsx", sheet = 2, skip = 23)

# load in metadata from Lupien lab
l_meta <- fread("metadata/brca.assignment.share.txt", nrows = 75) #rest of rows r empty


###########################################################
# Get sample names from tcga meta
###########################################################

meta <- meta[meta$case_id %in% l_meta$case_id,]
l_meta$sample_name <- meta[match(l_meta$case_id, meta$case_id),]$submitter_id
# note: b0f8d698-a30e-4d8d-b0a2-a5a01fac8406 | TCGA.A2.A0T4 is a duplicate

###########################################################
# Save metadata file
###########################################################

l_meta$sample_name <- gsub("-", "\\.", l_meta$sample_name)
write.csv(l_meta, file = "metadata/procdata/TCGA_subtype_label.csv", quote = F, row.names = F)


###########################################################
# Load in RNA data
###########################################################

# read in tcga tumour gene counts
counts <- read.table("data/rawdata/tcga/Human__TCGA_BRCA__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct", header = T)
counts <- counts[,colnames(counts) %in% l_meta$sample_name]
counts_mat <- t(counts)

# save
write.table(counts_mat, file = "data/procdata/TCGA/TCGA_BRCA_gene_counts.matrix", quote = F, sep = "\t", col.names = T)
