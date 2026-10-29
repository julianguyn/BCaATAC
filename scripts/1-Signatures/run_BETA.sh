#!/bin/bash

cd /Users/julianguyen/Documents/BETA2

###########################################################
# Run BETA for each ARCHE
###########################################################

for i in 1 2 3 4 5 6; do
    beta basic \
        -p ../BCaATAC/data/procdata/ARCHEs/BETAbeds/ARCHE${i}_50k.bed \
        -e ../BCaATAC/data/results/data/2-MolecularSigAnalysis/DEGforBETA/ARCHE${i}.txt \
        -k BSF \
        --gname2 \
        -g hg38 \
        --df 0.05 \
        -n ARCHE${i} \
        -o ../BCaATAC/data/results/data/1-Signatures/BETA
done
