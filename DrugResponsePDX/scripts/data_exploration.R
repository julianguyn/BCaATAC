

suppressMessages(library(Xeva))

# load in drug response data
#load("DrugResponsePDX/data/Xeva_PDXE.rds")
drug_screen <- read.csv("DrugResponsePDX/data/TNBC/drug_screening.csv")
# there is a column model.id

# from TNBC
models <- read.csv("DrugResponsePDX/data/TNBC/model_information.csv")
models2 <- read.csv("DrugResponsePDX/data/TNBC/batch_information.csv")

# from McGill
#models <- read.csv("DrugResponsePDX/data/McGill_TNBC/model_information.csv")

# load in Tina's samples
samples <- read.table("DrugResponsePDX/data/samplesTina/samples.tsv", header = T)
samples2 <- read.table("DrugResponsePDX/data/samplesTina/samples_organoids_run_bulkATAC.tsv", header = T)


# overlap
samples[samples$Sample %in% models$patient.id,]
# only 9 overlap


#################################################################################3
# read in samples we already have
samples <- read.csv("MetaData/Lupien/BCa_samples.csv")
samples <- samples[samples$type == "pdx",]

# remove low quality samples
samples <- samples[samples$qc == "pass",]

# save new dataframe
pdx_samples <- samples[,1:2]

# load in Tina's samples
samples <- read.table("DrugResponsePDX/data/samplesTina/samples.tsv", header = T)
tina <- samples$Sample
tina <- tina[order(tina)]

# check if samples are processed by Tina
pdx_samples$tina <- c(T,T,T,T,T,T,T,T,T,T,T,T,T,F,T,T,T,T,T,T,T,T,T,T,F,T,T)
pdx_samples$tina_name <- c("108099P1","16720","48602","53782","73720","BPTO95","BXTO143","BXTO152",
                           "DCBPTO19","DCBXTO147","DCBXTO137","DCBXTO138","DCBXTO66",NA,"INS_B014",
                           "INS_B019","REF_S_038","REF019","REF023","REF024","REF032","REF042",
                           "REF043","REF044",NA,"REFS034P1","REFS036P1A")


# check for overlap with TNBC drug response data
#pdx_samples$treatment <- c(T,F,T,T,F,F,F,F,F,F,F,F,F,F,T,T,F,T,T,T,F,F,F,F,F,F,T)
#pdx_samples$treatment_name <- c("108099", NA, "48602", "53782", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "INSB014","INSB019", NA, "REF019", "REF023", "REF024", NA, NA, NA, NA, NA, NA, "REFS036")


## read in all unique patient id's from XevaDB (from Matthew)
xeva <- read.csv("DrugResponsePDX/data/unique-patients.csv")

# check for overlap with Xeva samples
# check with Samah's list from Xeva: https://docs.google.com/spreadsheets/d/1neaM3TW01Ua1zVV1xN26X63YwW9IlGd5_UCdBJCfLfI/edit#gid=0
pdx_samples$xeva <- c("public","public","public",
                      "public","missing","public",
                      "private-OCTANE","missing","missing",
                      "missing","private-OCTANE","private-OCTANE",
                      "missing","missing","public",
                      "public","private","public",
                      "public","public","public",
                      "missing","private","public",
                      "public","public","public")
pdx_samples$xeva_name <- c("108099", "16720", "48602", 
                           "53782", NA, "BPTO.95", 
                           "BXTO.143", NA, NA, 
                           NA, "DCBXTO.137", "DCBXTO.138", 
                           NA, NA, "INSB014",
                           "INSB019", "REF-S-038", "REF019", 
                           "REF023", "REF024", "REFS032", 
                           NA, "REF043", "REF 044", 
                           "REF047", "REF-S-034", "REFS036")

# check if mapping for REFS032 (original: REF032_S18) is right
# DCBPTO19 might map to REF-S-034 DCPTOB19(WESONLY) in Samah's sheet

# add in Tina unique samples
tina_samples <- data.frame(sample = c(NA, NA, NA, NA),
                           filename = c(NA, NA, NA, NA),
                           tina = c(T,T,T,T),
                           tina_name = c("70420", "DCBXTO58", "BPTO93", "DCBXTO28"),
                           xeva = rep("missing", 4),
                           xeva_name = c(NA, NA, NA, NA))
pdx_samples <- rbind(pdx_samples, tina_samples)



write.table(pdx_samples, file = "DrugResponsePDX/data/pdx_available.tsv", quote = F, sep = "\t", row.names = F, col.names = T)

# get missing
x <- pdx_samples[pdx_samples$xeva == FALSE,]
write.table(x, file = "DrugResponsePDX/data/missing_from_Xeva.tsv", quote = F, sep = "\t", row.names = F, col.names = T)
