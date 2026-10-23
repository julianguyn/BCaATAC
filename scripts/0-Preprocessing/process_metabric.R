# load libraries
suppressPackageStartupMessages({
    library(data.table)
})

###########################################################
# Load in RNA data
###########################################################

# read in metabric illumina microarray expression (log2 intensity)
counts <- fread("data/rawdata/metabric/data_mrna_illumina_microarray.txt", data.table = FALSE)
counts$Entrez_Gene_Id <- NULL

###########################################################
# Collapse duplicated gene symbols
###########################################################

# get mean for hugo symbols that map to multiple probes
counts <- counts[!is.na(counts$Hugo_Symbol) & counts$Hugo_Symbol != "",]
counts <- as.data.table(counts)[, lapply(.SD, mean, na.rm = TRUE), by = Hugo_Symbol]
counts <- as.data.frame(counts)
counts[is.na(as.matrix(counts))] <- NA # mean of all NA returns NaN

###########################################################
# Format to match TCGA file
###########################################################

# order genes alphabetically to match tcga file
counts <- counts[order(counts$Hugo_Symbol),]

# rename gene column and sample names
colnames(counts)[1] <- "attrib_name"
colnames(counts) <- gsub("-", "\\.", colnames(counts))

# round values to 4 decimal places
counts[,-1] <- round(counts[,-1], 4)

# save
write.table(counts, file = "data/procdata/metabric/METABRIC_BRCA_Illumina_HT12v3_microarray_log2.cct",
            quote = F, sep = "\t", row.names = F, na = "NA")
