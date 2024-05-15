#!/usr/bin/env bash

set -euo pipefail

output_folder="resources"  # Modify this

wget \
    --continue \
    https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz


tar \
    -xvzf \
    gtdbtk_r214_data.tar.gz \
    -C $output_folder \
    --strip 1

rm -i gtdbtk_r214_data.tar.gz


echo "Remember to set $output_folder in config/features.yml"
