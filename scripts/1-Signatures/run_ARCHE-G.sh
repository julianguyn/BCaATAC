#!/bin/bash

cd /Users/julianguyen/Documents/BCaATAC

###########################################################
# Apply trained ARCHE-G surrogate models to new cohorts
###########################################################

OUT="data/results/data/1-Signatures/geneSurrogates"

# METABRIC
python3 scripts/1-Signatures/apply_ARCHE-G_surrogates.py \
    data/procdata/metabric/METABRIC_BRCA_Illumina_HT12v3_microarray_log2.cct \
    $OUT/metabric
