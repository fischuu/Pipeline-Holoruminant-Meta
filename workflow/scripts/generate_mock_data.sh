#!/usr/bin/env bash

set -euo pipefail

# Process chicken

mkdir --parents resources/reference/

wget \
    --continue \
    --output-document resources/reference/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.39.fa.gz \
    https://ftp.ensembl.org/pub/release-110/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.39.fa.gz

gzip -d resources/reference/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.39.fa.gz
bgzip -l9 -@ 8 resources/reference/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.39.fa

samtools faidx \
    resources/reference/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.39.fa.gz \
    39:100000-110000 \
| cut -f 1 -d ":" \
| bgzip -l 9 \
> resources/reference/chicken_39_sub.fa.gz

rm resources/reference/Gallus_gallus.* -f

# Get the mags and force to split them into 20kb contigs. binners and metabinners will complain for having a
# perfect simulation

seqtk seq resources/reference/mags.fa.gz -l 50000 \
| grep -v ^">" \
| awk '{print ">contig_" NR "\n" $1}' \
| pigz -11 > resources/reference/contigs.fa.gz


# Simulate mag reads
mkdir --parents resources/reads/

# Simulate one experiment. In the samples.tsv will trick it into thinking there are more
wgsim \
    -S 0 \
    -N 100000 \
     resources/reference/contigs.fa.gz \
    >(pigz -1 > resources/reads/sample1_1.fq.gz) \
    >(pigz -1 > resources/reads/sample1_2.fq.gz) \
> /dev/null

wgsim \
    -S 1 \
    -N 600000 \
    resources/reference/contigs.fa.gz \
    >(pigz -1 > resources/reads/sample2_1.fq.gz) \
    >(pigz -1 > resources/reads/sample2_2.fq.gz) \
> /dev/null

# 10M works
# 08M works
# 06M works
# 04M works
# 02M works
# 01M works
# 500k does not work
# 600k works
# 550k does not work
