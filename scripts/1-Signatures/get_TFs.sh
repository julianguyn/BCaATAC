#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH -c 1
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=jlm.nguyen@mail.utoronto.ca
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --job-name=subset_tf
#SBATCH --output=logs/subset_tf.out

# cd projects/BCaATAC/Griffin/data/procdata/TF

INTERVAL_FILE="Homo_sapiens_meta_clusters.interval"
GENE_LIST="homer_TF_genes.txt"
OUTPUT_DIR="ARCHE_TFs"

mkdir -p "$OUTPUT_DIR"

while IFS= read -r gene || [ -n "$gene" ]; do
    [ -z "$gene" ] && continue
    echo "$gene"
    OUTFILE="${OUTPUT_DIR}/${gene}.bed"
    printf "seqnames\tstart\tend\n" > "$OUTFILE"
    awk -v gene="$gene" 'BEGIN {OFS=FS="\t"} $6 == gene {chrom=$1; sub(/^chr/, "", chrom); print chrom, $2, $3}' "$INTERVAL_FILE" >> "$OUTFILE"
done < "$GENE_LIST"