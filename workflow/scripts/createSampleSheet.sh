#!/bin/bash

# Directory containing the .fastq.gz files
fastq_path="reads/"

# Output file for the samplesheet
output="samplesheet.tsv"

# Adapters
forward_adapter="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
reverse_adapter="CTGTCTCTTATACACATCTGACGCTGCCGACGA"

# Initialize the samplesheet with headers
echo -e "sample_id\tlibrary_id\tforward_filename\treverse_filename\tforward_adapter\treverse_adapter\tassembly_ids" > $output

# Find all .fastq.gz files and process them
for file in *_R1_001.fastq.gz; do
    if [[ -f "$file" ]]; then
        sample_id=$(basename "$file" | sed -E 's/(.*)_S[0-9]+_L[0-9]{3}_R1_001.fastq.gz/\1/')
        forward_file="${fastq_path}$(basename "$file")"
        reverse_file="${fastq_path}$(basename "${file/_R1_001.fastq.gz/_R2_001.fastq.gz}")"
        
        # Check if the reverse file exists
        if [[ -f "${file/_R1_001.fastq.gz/_R2_001.fastq.gz}" ]]; then
            echo -e "$sample_id\tlib1\t$forward_file\t$reverse_file\t$forward_adapter\t$reverse_adapter\t$sample_id" >> $output
        fi
    fi
done

echo "Samplesheet created: $output"

