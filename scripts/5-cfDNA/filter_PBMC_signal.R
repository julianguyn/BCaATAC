# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(ggplot2)
})

source("utils/palettes.R")
source("utils/filter_PBMC_signal.R")

###########################################################
# Load in differential peaks
###########################################################

dir <- "data/rawdata/cfDNA/Collaboration_w_Tina_Griffin_Reference/"

# blood enriched shared overlap any
filename <- "Blood_enriched_shared_overlap_any.bed"
blood_all <- read.table(paste0(dir, filename))
colnames(blood_all) <- c("chr", "start", "end")

#blood enriched vs ER
filename <- "Blood_enriched_vs_ER_FDR005_FOLD1.bed"
blood_er <- read.table(paste0(dir, filename))
colnames(blood_er) <- c("chr", "start", "end")

###########################################################
# Find hits in ARCHE beds
###########################################################

filter_PBMC_signal(blood_er, "BloodvsER")
filter_PBMC_signal(blood_all, "BloodvsAll")

###########################################################
# TODO: compare sites filtered across different differntials
###########################################################