# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(ggplot2)
    library(ComplexHeatmap)
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

# blood enriched vs Basal
filename <- "Blood_enriched_vs_Basal_FDR005_FOLD1.bed"
blood_basal <- read.table(paste0(dir, filename))
colnames(blood_basal) <- c("chr", "start", "end")

#blood enriched vs HER2
filename <- "Blood_enriched_vs_Her2_FDR005_FOLD1.bed"
blood_her2 <- read.table(paste0(dir, filename))
colnames(blood_her2) <- c("chr", "start", "end")

###########################################################
# Find hits in ARCHE beds
###########################################################

removed_all <- filter_PBMC_signal(blood_all, "BloodvsAll")
removed_er <- filter_PBMC_signal(blood_er, "BloodvsER")
removed_basal <- filter_PBMC_signal(blood_basal, "BloodvsBasal")
removed_her2 <- filter_PBMC_signal(blood_her2, "BloodvsHer2")

###########################################################
# Compare sites filtered across different cases
###########################################################

plot_sites_UPSET("ARCHE1", arche_list = TRUE)
plot_sites_UPSET("ARCHE2", arche_list = TRUE)
plot_sites_UPSET("ARCHE3", arche_list = TRUE)
plot_sites_UPSET("ARCHE4", arche_list = TRUE)
plot_sites_UPSET("ARCHE5", arche_list = TRUE)
plot_sites_UPSET("ARCHE6", arche_list = TRUE)

plot_sites_UPSET(removed_all, "All_Subtypes")
plot_sites_UPSET(removed_er, "ER")
plot_sites_UPSET(removed_basal, "Basal")
plot_sites_UPSET(removed_her2, "Her2")
