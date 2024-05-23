# Script to extract gene counts matrix from UHN Breast PSet

setwd("C:/Users/julia/Documents/BCaATAC")

suppressMessages(library(PharmacoGx))

# load in cell line RNA seq counts
uhnbreast2 <- readRDS("DrugResponse/data/PharmacoSet.RDS")
cells <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@assays@data$expr
colnames(cells) <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@colData$sampleid
genes_meta <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@rowRanges
rownames(cells) <- uhnbreast2@molecularProfiles@ExperimentList$genes_counts@rowRanges$gene_name

# load in cell line annotation
cells_meta <- read.csv("MetaData/Lupien/BCa_samples.csv")
# from https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2911-z/tables/3, set 600MPE to LuminalA
cells_meta[which(cells_meta$subtype == "cell_line"),]$subtype <- "LumA"

# keep only cell lines being used
samples <- read.csv("DrugResponse/data/cl_samples.csv")
samples$file <- gsub("\\.$", "", samples$file)
samples$file <- gsub("\\.", "-", samples$file)
samples$match <- gsub("-", "", samples$sample)

# remove dups
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz

# save subtype
cells_meta <- cells_meta[which(cells_meta$sample %in% samples$file),]
samples$subtype <- cells_meta[match(samples$file, cells_meta$sample),]$subtype

# match cells to metadata (samples)
# missing: "SUM52PE" "HBL100" "HCC1008" 
colnames(cells)[which(colnames(cells) == "LY2")] <- "MCF7/LY2"
colnames(cells)[which(colnames(cells) == "SUM149")] <- "SUM149PT"
colnames(cells)[which(colnames(cells) == "SUM159")] <- "SUM159PT"
colnames(cells)[which(colnames(cells) == "HS578T")] <- "Hs 578T"

# save only needed cells
cells <- cells[,colnames(cells) %in% samples$match]

# format dataframe
cells <- as.data.frame(t(cells))

write.table(cells, file = "Signatures/data/bcacells_gene_counts.matrix",  quote = F, sep = "\t", col.names = T)
