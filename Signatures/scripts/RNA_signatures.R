setwd("C:/Users/julia/Documents/BCaATAC")


library(data.table)
library(readxl)

# load in annotation file from https://linkedomics.org/data_download/TCGA-BRCA/
tcga <- read.table("Signatures/data/Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi", header = T)
tcga <- colnames(tcga)[2:1098]


# map sample file case_id to Case ID from https://portal.gdc.cancer.gov/v1/annotations
# found annotation spreadsheet
meta <- read_excel("Signatures/data/TCGA-ATAC_DataS1_DonorsAndStats_v4.xlsx", sheet = 2, skip = 23)

lupien <- fread("Signatures/data/brca.assignment.share.txt")
lupien <- lupien[1:75,]
meta <- meta[meta$case_id %in% lupien$case_id,]
lupien$sample_name <- meta[match(lupien$case_id, meta$case_id),]$submitter_id
# note: b0f8d698-a30e-4d8d-b0a2-a5a01fac8406 | TCGA.A2.A0T4 is a duplicate

lupien$sample_name <- gsub("-", "\\.", lupien$sample_name)
write.csv(lupien, file = "Signatures/data/TCGA_subtype_label.csv", quote = F, row.names = F)

# read in tcga tumour gene counts from https://linkedomics.org/data_download/TCGA-BRCA/
counts <- read.table("Signatures/data/Human__TCGA_BRCA__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct", header = T)
counts <- counts[,colnames(counts) %in% lupien$sample_name]

counts_mat <- t(counts)
write.table(counts_mat, file = "Signatures/data/TCGA_BRCA_gene_counts.matrix", quote = F, sep = "\t", col.names = T)
