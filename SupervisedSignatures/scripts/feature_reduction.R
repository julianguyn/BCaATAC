# Identify peaks significantly associated with drug response (univariable)

setwd("C:/Users/julia/Documents/BCaATAC")

library(data.table)
suppressMessages(library(readxl))
suppressMessages(library(survcomp))

### Identify drugs of interest ###

# read in PDX drug response data
pdx_response <- as.data.frame(read_excel("DrugResponsePDX/data/drugresponse/DrugResponse_PDX.xlsx", sheet = 1))
pdx_response[pdx_response == "NA"] <- NA
pdx_response[pdx_response == "TAXOL"] <- "PACLITAXEL"
pdx_response[pdx_response == "CARBOPLATIN"] <- "CARBOPLATINUM"
pdx_response <- na.omit(pdx_response)

# remove CONTROLs and H2O samples (left with n=23 drugs)
drugs_of_interest <- unique(pdx_response$drug[-which(pdx_response$drug %in% c("H2O", pdx_response$drug[grep("CONTROL", pdx_response$drug)]))])
pdx_response <- pdx_response[pdx_response$drug %in% drugs_of_interest,]

# keep only drugs with at least 10 samples (left with n=11 drugs)
drugs_of_interest <- names(table(pdx_response$drug)[table(pdx_response$drug) > 10])
pdx_response <- pdx_response[pdx_response$drug %in% drugs_of_interest,]


### Get cell line drug response for drugs of interest ### 

# load in cell line drug response data
load("DrugResponse/results/data/sensitivity_data.RData")

# load in and extract needed cell lines
samples <- read.csv("DrugResponse/data/cl_samples.csv")
dup <- samples$sample[duplicated(samples$sample)]
tmp <- samples[samples$sample %in% dup,]
samples <- rbind(samples[-which(samples$sample %in% dup),], tmp[which(tmp$seq == "Nergiz"),]) # from dups, keep only the ones by Nergiz


# keep only drugs of interest
ubr1_sen <- ubr1_sen[toupper(rownames(ubr1_sen)) %in% drugs_of_interest,colnames(ubr1_sen) %in% samples$sample] #taxol only, 43 cells
ubr2_sen <- ubr2_sen[toupper(rownames(ubr2_sen)) %in% drugs_of_interest,colnames(ubr2_sen) %in% samples$sample] #taxol, carboplatin, & eribulin, 42 cells
gray_sen <- gray_sen[toupper(rownames(gray_sen)) %in% drugs_of_interest,colnames(gray_sen) %in% samples$sample] #taxol & carboplatin, 41 cells
gcsi_sen <- gcsi_sen[toupper(rownames(gcsi_sen)) %in% drugs_of_interest,colnames(gcsi_sen) %in% samples$sample] #taxol only, 25 cells
gdsc_sen <- gdsc_sen[toupper(rownames(gdsc_sen)) %in% drugs_of_interest,colnames(gdsc_sen) %in% samples$sample] #taxol only, 34 cells
ctrp_sen <- ctrp_sen[toupper(rownames(ctrp_sen)) %in% drugs_of_interest,colnames(ctrp_sen) %in% samples$sample] #taxol & carboplatin, 32 cells
ccle_sen <- ccle_sen[toupper(rownames(ccle_sen)) %in% drugs_of_interest,colnames(ccle_sen) %in% samples$sample] #taxol only, 23 cells


# save to H4H 
# save(samples, ubr1_sen, ubr2_sen, gray_sen, gcsi_sen, gdsc_sen, ctrp_sen, ccle_sen, file = "temp_toH4H.RData")

######## TEMPORARY, FOR H4H
library(data.table)
suppressMessages(library(survcomp))

load("../temp_toH4H.RData")
##############################

# colnames of sensitivity has to match with peaks

# load in peak calls
peak_calls <- as.data.frame(fread("../bcacell_lines.tsv"))
rownames(peak_calls) <- paste(peak_calls$V1, peak_calls$V2, peak_calls$V3, sep = ":")
peak_calls <- peak_calls[,-c(1:3)]

# keep only needed samples and rename sample names
colnames(peak_calls) <- gsub("-", "\\.", colnames(peak_calls))
peak_calls <- peak_calls[,which(colnames(peak_calls) %in% samples$file)]
colnames(peak_calls) <- samples[match(colnames(peak_calls), samples$file),]$sample


### Compute concordance index to find univariable associations ###

# function to compute concordance index
computeCI <- function(signature_scores, sensitivity_data, label) {

    # create data frame to hold results
    combinations <- as.data.frame(matrix(data = NA, nrow = nrow(signature_scores) * nrow(sensitivity_data), ncol = 7))
    colnames(combinations) <- c("peak", "drug", "ci", "pvalue", "se", "upper", "lower")
    combinations$peak <- rep(rownames(signature_scores), nrow(sensitivity_data))
    combinations$drug <- rep(rownames(sensitivity_data), each = nrow(signature_scores))

    # compute concordance index
    for (i in 1:nrow(combinations)){
        if (i %% 10000 == 0) {print(paste0(i, " out of ", nrow(combinations), " complete"))}

        ci <- survcomp::concordance.index(as.numeric(sensitivity_data[combinations[,2][i],]), 
                                            # row from sensitivity_data with the drug for all samples
                                            surv.time = as.numeric(unlist(-signature_scores[combinations[,1][i],])), 
                                            # row of signature i in signature scores matrix with cell lines as columns
                                            surv.event = rep(1,length(sensitivity_data)), 
                                            # df of all drugs as rows with all samples as columns
                                            outx = TRUE, method="noether", na.rm = TRUE)

        combinations$pvalue[i] <- ci$p.value
        combinations$ci[i] <- ci$c.index
        combinations$se[i] <- ci$se
        combinations$upper[i] <- ci$upper
        combinations$lower[i] <- ci$lower

        
    }

    # filtering and multiple test correction
    combinations <- combinations[complete.cases(combinations$pvalue),]
    combinations$FDR <- p.adjust(combinations$pvalue, method = "BH", n = length(combinations$pvalue))
    combinations$FDRsig <- ifelse(combinations$FDR < 0.05, TRUE, FALSE)

    
    # format dataframe for plotting (order by ci and add variable rank)
    combinations <- combinations[order(combinations$ci),]
    combinations$rank <- 1:nrow(combinations)
    combinations$pairs <- paste0(combinations$peak, "-", combinations$drug)
    combinations$pset <- c(rep(label, nrow(combinations)))
    
    return(combinations)
}

# subset peak calls to samples in dataset
ubr1_sig <- peak_calls[,which(colnames(peak_calls) %in% colnames(ubr1_sen))]
ubr2_sig <- peak_calls[,which(colnames(peak_calls) %in% colnames(ubr2_sen))]
gray_sig <- peak_calls[,which(colnames(peak_calls) %in% colnames(gray_sen))]
gcsi_sig <- peak_calls[,which(colnames(peak_calls) %in% colnames(gcsi_sen))]
gdsc_sig <- peak_calls[,which(colnames(peak_calls) %in% colnames(gdsc_sen))]
ctrp_sig <- peak_calls[,which(colnames(peak_calls) %in% colnames(ctrp_sen))]
ccle_sig <- peak_calls[,which(colnames(peak_calls) %in% colnames(ccle_sen))]

# compute associations
ubr1_com <- computeCI(ubr1_sig, ubr1_sen, "uhnbreast1")
ubr2_com <- computeCI(ubr2_sig, ubr2_sen, "uhnbreast2")
gray_com <- computeCI(gray_sig, gray_sen, "gray")
gcsi_com <- computeCI(gcsi_sig, gcsi_sen, "gcsi")
gdsc_com <- computeCI(gdsc_sig, gdsc_sen, "gdsc")
ctrp_com <- computeCI(ctrp_sig, ctrp_sen, "ctrp")
ccle_com <- computeCI(ccle_sig, ccle_sen, "ccle")



save(ubr1_com, file = "../results/res1.RData")
save(ubr2_com, file = "../results/res2.RData")
save(gray_com, file = "../results/res3.RData")
save(gcsi_com, file = "../results/res4.RData")
save(gdsc_com, file = "../results/res5.RData")
save(ctrp_com, file = "../results/res6.RData")
save(ccle_com, file = "../results/res7.RData")

# save(ubr1_com, ubr2_com, gray_com, gcsi_com, gdsc_com, ctrp_com, ccle_com, file = "SupervisedDrugResponseAssociations.RData")


