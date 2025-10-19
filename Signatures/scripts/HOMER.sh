#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=4:00:00

#SBATCH --job-name=runHOMER
#SBATCH --output=logs/runHOMER.out


###########################################################
# Specify genome and annotation file
###########################################################

GENOME="references/GENCODE/GRCh38_v45/genome.fa"
GTF="references/GENCODE/GRCh38_v45/annotation.gtf"

###########################################################
# Annotate ARCHE peak sets
###########################################################

annotatePeaks.pl ../data/ARCHE1_20k.bed $GENOME -gtf $GTF -size given > ../results/annotatePeaks/ARCHE1.txt
annotatePeaks.pl ../data/ARCHE2_20k.bed $GENOME -gtf $GTF -size given > ../results/annotatePeaks/ARCHE2.txt
annotatePeaks.pl ../data/ARCHE3_20k.bed $GENOME -gtf $GTF -size given > ../results/annotatePeaks/ARCHE3.txt
annotatePeaks.pl ../data/ARCHE4_20k.bed $GENOME -gtf $GTF -size given > ../results/annotatePeaks/ARCHE4.txt
annotatePeaks.pl ../data/ARCHE5_20k.bed $GENOME -gtf $GTF -size given > ../results/annotatePeaks/ARCHE5.txt
annotatePeaks.pl ../data/ARCHE6_20k.bed $GENOME -gtf $GTF -size given > ../results/annotatePeaks/ARCHE6.txt

###########################################################
# Motif enrichment
###########################################################

findMotifsGenome.pl ../data/ARCHE1_20k.bed $GENOME -gtf $GTF ../results/findMotifsGenome/ARCHE1 -bg ../data/Background.bed -size given
findMotifsGenome.pl ../data/ARCHE2_20k.bed $GENOME -gtf $GTF ../results/findMotifsGenome/ARCHE2 -bg ../data/Background.bed -size given
findMotifsGenome.pl ../data/ARCHE3_20k.bed $GENOME -gtf $GTF ../results/findMotifsGenome/ARCHE3 -bg ../data/Background.bed -size given
findMotifsGenome.pl ../data/ARCHE4_20k.bed $GENOME -gtf $GTF ../results/findMotifsGenome/ARCHE4 -bg ../data/Background.bed -size given
findMotifsGenome.pl ../data/ARCHE5_20k.bed $GENOME -gtf $GTF ../results/findMotifsGenome/ARCHE5 -bg ../data/Background.bed -size given
findMotifsGenome.pl ../data/ARCHE6_20k.bed $GENOME -gtf $GTF ../results/findMotifsGenome/ARCHE6 -bg ../data/Background.bed -size given





# modifications to scripts to above two and preparseGenome.pl and assignGeneWeights.pl and compareMotifs.pl
use lib "/cluster/projects/bhklab/projects/BCaATAC/HOMER/.//bin";
my $homeDir = "/cluster/projects/bhklab/projects/BCaATAC/HOMER/./";

my $genomeDir = "/cluster/projects/bhklab/projects/BCaATAC/HOMER/data/genomes/hg38";
my $preparsedDir = "/cluster/projects/bhklab/projects/BCaATAC/HOMER/data/genomes/hg38/preparsed/";
my $genomeParseDir = "/cluster/projects/bhklab/projects/BCaATAC/HOMER/data/genomes/hg38/preparsed/";
my $preparsedDirFromConfig = "/cluster/projects/bhklab/projects/BCaATAC/HOMER/data/genomes/hg38/preparsed/";
my $defFile = "/cluster/projects/bhklab/projects/BCaATAC/HOMER/data/genomes/hg38/hg38.tss";

check:
findKnownMotifs.pl
compareMotifs.pl
makeBinaryFile.pl
adjustRedunGroupFile.pl
randRemoveBackground.pl

