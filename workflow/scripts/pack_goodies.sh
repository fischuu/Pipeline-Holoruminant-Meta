#!/usr/bin/env bash

set -euo pipefail

tar \
    --create \
    --verbose \
    --file final_results.tar.gz \
    --compress-program="pigz" \
    reports \
    results/dereplicate/drep/dereplicated_genomes.fa \
    results/dereplicate/coverm/*.tsv \
    results/dereplicate/checkm2/quality_report.tsv \
    results/dereplicate/gtdbtk/gtdbtk.summary.tsv \
    results/dereplicate/dram/annotations.tsv \
    results/dereplicate/dram/product.tsv \
    results/dereplciate/dram/product.html \
    results/dereplicate/dram/metabolism_summary.xlsx
