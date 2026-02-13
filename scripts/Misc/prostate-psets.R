# prostate cancer cell lines (for Lisanne)

# load libraries
suppressPackageStartupMessages({
  library(PharmacoGx)
})

###########################################################
# Download PSets
###########################################################

# to see a list of other available PSets, you can use
# `availablePSets()`

gdsc <- downloadPSet("GDSC_2020(v2-8.2)")
ccle <- downloadPSet("CCLE_2015")
ctrp <- downloadPSet("CTRPv2_2015")
gcsi <- downloadPSet("gCSI_2019")

###########################################################
# Check number of prostate cell lines
###########################################################

# GDSC: 9 prostate cell lines
table(gdsc@sample$tissueid) 
rownames(gdsc@sample[gdsc@sample$tissueid == "Prostate",])
# BPH-1, DU145, LNCaP clone FGC, PC-3, 22Rv1, VCaP, NCI-H660, RWPE-1, RWPE2-W99

# CCLE: 8 prostate cell lines
table(ccle@sample$tissueid)
rownames(ccle@sample[ccle@sample$tissueid == "Prostate",])
# 22Rv1, DU145, LNCaP clone FGC, MDA-PCa-2b, NCI-H660, PC-3, PrEC LH, VCaP

# gCSI: 5 prostate cell lines
table(gcsi@sample$tissueid)
rownames(gcsi@sample[gcsi@sample$tissueid == "Prostate",])
# LNCaP clone FGC, DU145, PC-3, LNCaP, 22Rv1

###########################################################
# Get RNA
###########################################################

# this gives you a summarized experiment with metadata
rna <- summarizeMolecularProfiles(
    gdsc, # replace with pset name
    mDataType = "Kallisto_0.46.1.rnaseq.counts"
)

# this gives you just the rna counts matrix
rna_counts <- assay(rna)

###########################################################
# Get drug response
###########################################################

# this gives you a drug response matrix (values are AAC)
sensitivity <- summarizeSensitivityProfiles(
    gdsc, # replace with pset name
    sensitivity.measure = "aac_recomputed",
    summary.stat = "median"
)
