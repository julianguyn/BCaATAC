#!/bin/bash

cd /Users/julianguyen/Documents/BCaATAC

###########################################################
# Apply trained ARCHE-G surrogate models to new cohorts
###########################################################

OUT="data/results/data/1-Signatures/geneSurrogates"

# PDXE
python3 scripts/1-Signatures/GenomicSurrogates/apply_ARCHE-G_surrogates.py \
    data/procdata/PDXs/pdxe/rna.csv \
    rna_model \
    $OUT/pdxe

# METABRIC
python3 scripts/1-Signatures/GenomicSurrogates/apply_ARCHE-G_surrogates.py \
    data/procdata/metabric/METABRIC_BRCA_Illumina_HT12v3_microarray_log2.cct \
    microarray_model \
    $OUT/metabric
