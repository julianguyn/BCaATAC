# load libraries
suppressPackageStartupMessages({
  library(reshape2)
  library(readxl)
  library(ggplot2)
  library(data.table)
  library(dplyr)
})

# ---------------------------------------------------------
# Make samples.yaml file
# ---------------------------------------------------------

# compiled from H4H
load("data/temp/bams.RData")

# format
df <- data.frame(id = NA, path = bams)
df$id <- gsub(".*/", "", df$path)
df$id <- gsub("\\..*", "", df$id)
df$id <- paste0("  ", df$id, ":")

write.table(df, file = "data/temp/reflect_samples.txt", quote = FALSE, row.names = FALSE)


###########################################################
# Load in data
###########################################################

# load in metadata
meta <- read_excel("data/rawdata/cfdna/REFLECT-6B_metadata.xlsx", sheet = 1) |> as.data.frame()
