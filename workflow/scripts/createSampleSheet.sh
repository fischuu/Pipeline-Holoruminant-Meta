#!/bin/bash

# Directory containing the .fastq.gz files
fastq_path="reads/"

# Output file for the samplesheet
output="config/samples.tsv"

# Adapters
forward_adapter="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
reverse_adapter="CTGTCTCTTATACACATCTGACGCTGCCGACGA"

# Initialize the samplesheet with headers
echo -e "sample_id\tlibrary_id\tforward_filename\treverse_filename\tforward_adapter\treverse_adapter\tassembly_ids" > $output

echo "Looking for files in ${fastq_path}..."

# Find all _R1_ fastq.gz files and process them
found_files=0
for file in ${fastq_path}*_R1_*.fastq.gz; do
    if [[ -f "$file" ]]; then
        found_files=1
        # Extract sample ID (everything before the first underscore)
        sample_id=$(basename "$file" | cut -d'_' -f1)
        forward_file="${fastq_path}$(basename "$file")"
        reverse_file="${fastq_path}$(basename "${file/_R1_/_R2_}")"
        
        # Check if the reverse file exists
        if [[ -f "$reverse_file" ]]; then
            echo "Processing sample: $sample_id"
            echo -e "$sample_id\tlib1\t$forward_file\t$reverse_file\t$forward_adapter\t$reverse_adapter\t$sample_id" >> $output
        else
            echo "Reverse file not found for sample: $sample_id (Expected: $reverse_file)"
        fi
    else
        echo "No forward files found matching pattern *_R1_*.fastq.gz"
    fi
done

if [[ $found_files -eq 0 ]]; then
    echo "No files found in ${fastq_path} matching the pattern *_R1_*.fastq.gz"
fi

echo "Samplesheet created: $output"
