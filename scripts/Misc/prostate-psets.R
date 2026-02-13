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
# Get drug response
###########################################################

# this gives you a drug response matrix (values are AAC)
gdsc_sen <- summarizeSensitivityProfiles(
    gdsc, # replace with pset name
    sensitivity.measure = "aac_recomputed",
    summary.stat = "median"
)
ctrp_sen <- summarizeSensitivityProfiles(
    ctrp, # replace with pset name
    sensitivity.measure = "aac_recomputed",
    summary.stat = "median"
)
gcsi_sen <- summarizeSensitivityProfiles(
    gcsi, # replace with pset name
    sensitivity.measure = "aac_recomputed",
    summary.stat = "median"
)

save(gdsc_sen, ctrp_sen, gcsi_sen, file = "temp/toLisanne.RData")

###########################################################
# Check for drugs in drug names
###########################################################

drugs_of_interest <- c(
    "Tazemetostat",
    "Valemetostat",
    "Tulmimetostat",
    "CPI-0209",
    "Igermetostat",
    "XNW5004",
    "TR115",
    "Mevrometostat",
    "GSK126",
    "EED226",
    "UNC1999",
    "Azacitidine",
    "Decitabine",
    "Cedazuridine",
    "Guadecitabine",
    "CC-486",
    "Zebularine",
    "GSK3685032"
)

drugs_of_interest[drugs_of_interest %in% rownames(ctrp_sen)]
drugs_of_interest[drugs_of_interest %in% rownames(gcsi_sen)]

ctrp_sen[
    drugs_of_interest[drugs_of_interest %in% rownames(ctrp_sen)],
    which(colnames(ctrp_sen) %in% rownames(ccle@sample[ccle@sample$tissueid == "Prostate",]))
]

gcsi_sen[
    drugs_of_interest[drugs_of_interest %in% rownames(gcsi_sen)],
    which(colnames(gcsi_sen) %in% rownames(gcsi@sample[gcsi@sample$tissueid == "Prostate",])),
    drop = FALSE
]