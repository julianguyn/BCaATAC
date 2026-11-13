# load libraries
suppressPackageStartupMessages({
  library(readxl)
  library(reshape2)
  library(stringr)
  library(Xeva)
  library(data.table)
  library(org.Hs.eg.db)
})

# ---------------------------------------------------------
# Process internal xeva data

###########################################################
# Load in data
###########################################################

auc <- read_excel("data/rawdata/pdx/auc.xlsx", sheet = 1)
rec <- read_excel("data/rawdata/pdx/mRECIST.xlsx", sheet = 1)

###########################################################
# Parse each replicate into one sample value
###########################################################

# helper function to melt and parse replicates
parse_reps <- function(df) {

    colnames(df)[1] <- 'drug'
    df <- reshape2::melt(df, id = 'drug')
    df <- df[!is.na(df$value),]

    # track rows to remove
    to_rm <- c()

    # iterate through each row, split replicates
    for (i in 1:nrow(df)) {

        split <-  str_split_1(df$value[i], ";")
        n_reps <- length(split)

        # if replicates exist
        if (n_reps > 1) {
            drug <- df$drug[i]
            var <- df$variable[i]

            # add each replicate as it's individual sample
            for (j in 1:n_reps) {
                df <- rbind(df, data.frame(drug = drug, variable = var, value = split[j]))
            }
            
            # track rows with replicates
            to_rm <- c(to_rm, i)
        }
    }
    
    # remove original rows with replicates
    df <- df[-c(to_rm),]
}

auc <- parse_reps(auc)
rec <- parse_reps(rec)


###########################################################
# Format dataframes
###########################################################

colnames(auc) <- c("drug", "patient.id", "AUC")
colnames(rec) <- c("drug", "patient.id", "mRECIST")

###########################################################
# Save files
###########################################################

write.csv(auc, file = "data/procdata/PDXs/drugresponse/auc_reps.csv", quote = F, row.names = F)
write.csv(rec, file = "data/procdata/PDXs/drugresponse/mRECIST_reps.csv", quote = F, row.names = F)


# ---------------------------------------------------------
# Process PDXE

pdxe <- readRDS("data/rawdata/pdx/Xeva_PDXE.rds")

###########################################################
# Process drug response
###########################################################

# get drug response
drug <- pdxe@model
bin <- drug[grep("binimetinib", drug$drug),]
model.ids <- bin$model.id[bin$oncotree.code == "BREAST"]
#tra <- drug[grep("trametinib", drug$drug),] # this is only PDAC and CRC
sen <- pdxe@sensitivity$model
sen <- sen[sen$model.id %in% model.ids,]
sen$model.id <- sub("\\.", "-", sub(".biib", "", sen$model.id))

###########################################################
# Process molecular profiles
###########################################################

# get molecular profiles
rna <- pdxe@molecularProfiles$RNASeq@assayData$exprs
mut <- pdxe@molecularProfiles$mutation@assayData$exprs

# keep samples (n=38)
keep <- intersect(colnames(rna), sen$model.id)
rna <- rna[,keep]
mut <- mut[,keep]
sen <- sen[match(keep, sen$model.id),]

###########################################################
# Match gene names to TCGA (for ARCHE-G surrogate models)
###########################################################

# rename PDXE genes to their TCGA symbol 
tcga_genes <- fread("data/rawdata/tcga/Human__TCGA_BRCA__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct", select = 1)[[1]]
old_names <- setdiff(tcga_genes, rownames(rna))

alias <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = old_names, keytype = "ALIAS", columns = "SYMBOL"))
alias <- alias[!is.na(alias$SYMBOL) & alias$SYMBOL %in% rownames(rna) & !alias$SYMBOL %in% tcga_genes,]
alias <- alias[!alias$ALIAS %in% alias$ALIAS[duplicated(alias$ALIAS)] &
               !alias$SYMBOL %in% alias$SYMBOL[duplicated(alias$SYMBOL)],] # one-to-one only

rownames(rna)[match(alias$SYMBOL, rownames(rna))] <- alias$ALIAS
message("Renamed ", nrow(alias), " PDXE genes to TCGA symbols")

###########################################################
# Save profiles
###########################################################

write.csv(rna, file = "data/procdata/PDXs/pdxe/rna.csv", quote = F, row.names = T)
write.csv(mut, file = "data/procdata/PDXs/pdxe/mut.csv", quote = F, row.names = T)
write.csv(sen, file = "data/procdata/PDXs/pdxe/sen.csv", quote = F, row.names = T)
