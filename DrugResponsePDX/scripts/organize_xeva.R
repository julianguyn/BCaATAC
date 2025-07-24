# script to organize treatment response files for xeva

setwd("C:/Users/julia/Documents/BCaATAC")

# load libraries
suppressPackageStartupMessages({
  library(readxl)
  library(reshape2)
  library(stringr)
})

###########################################################
# Load in data
###########################################################

auc <- read_excel("DrugResponsePDX/data/Jul232025-Xeva/auc.xlsx", sheet = 1)
rec <- read_excel("DrugResponsePDX/data/Jul232025-Xeva/mRECIST.xlsx", sheet = 1)


###########################################################
# Parse each replicate into one sample value
###########################################################

# define function to melt and parse replicates
parse_reps <- function(df) {

    colnames(df)[1] <- 'drug'
    df <- melt(df, id = 'drug')
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

write.csv(auc, file = "DrugResponsePDX/data/Jul232025-Xeva/auc_reps.csv", quote = F, row.names = F)
write.csv(rec, file = "DrugResponsePDX/data/Jul232025-Xeva/mRECIST_reps.csv", quote = F, row.names = F)