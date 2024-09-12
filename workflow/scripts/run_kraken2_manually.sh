#!/usr/bin/env bash

set -euo pipefail

db_path="resources/kraken2_mock"
db_name=$(basename $db_path)
shm_dbname="/dev/shm/$db_name"

fastp_dir="results/preprocess/fastp/"
kraken2_dir="kraken2"

mapfile -t sample_ids < <(find "$fastp_dir" -name "*_1.fq.gz" -exec basename {} _1.fq.gz \;)

rsync -Pravt $db_path /dev/shm/

mkdir --parents "$kraken2_dir"

for sample_id in "${sample_ids[@]}" ; do

    echo "Processing sample ${sample_id} ..."

    kraken2 \
        --db "$shm_dbname" \
        --threads 24 \
        --paired \
        --gzip-compressed \
        --output >(pigz > "$kraken2_dir/$sample_id.kraken2") \
        --report "$kraken2_dir/$sample_id.kraken2.report" \
        "$fastp_dir/${sample_id}_1.fq.gz" \
        "$fastp_dir/${sample_id}_2.fq.gz" \
    2> "$kraken2_dir/${sample_id}.kraken2.log" 1>&2

done

rm -rf "$shm_dbname"
