#!/usr/bin/env bash
set -euo pipefail

# Assume kraken2 is installed

# Build db


# Has to be a fasta
gzip --decompress --keep resources/reference/mags.fa.gz


kraken2-build \
    --no-masking \
    --add-to-library resources/reference/mags.fa \
    --db resources/kraken2_mock \
    --threads 4

rm -f  resources/reference/mags.fa

kraken2-build \
    --download-taxonomy \
    --db resources/kraken2_mock

kraken2-build \
    --build \
    --db resources/kraken2_mock

kraken2-build \
    --clean \
    --db resources/kraken2_mock

# Test

kraken2 \
    --db resources/kraken2_mock \
    --threads 4 \
    --report resources/kraken2_mock.report \
    --output resources/kraken2_mock.out \
    --gzip-compressed \
    resources/reference/mags.fa.gz


# Run

kraken2 \
    --db resources/kraken2_mock \
    --threads 8 \
    --report /dev/null \
    --output /dev/null \
    --gzip-compressed \
    --paired \
    resources/reads/sample1_1.fq.gz \
    resources/reads/sample1_2.fq.gz
