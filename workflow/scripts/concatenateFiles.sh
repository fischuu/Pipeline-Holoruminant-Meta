#!/bin/bash

# Check for dry-run option
dry_run=false
if [[ "$1" == "--dry-run" ]]; then
    dry_run=true
fi

# List all files in the current directory
files=$(ls *.fastq.gz)

# Loop through each unique sample prefix, adjust the cut for your data
for sample in $(echo "$files" | cut -d '_' -f 1 | sort | uniq); do

    # Find R1 and R2 files (L003 and L004)
    R1_files=$(echo "$files" | grep "^$sample" | grep "_R1_" | tr '\n' ' ')
    R2_files=$(echo "$files" | grep "^$sample" | grep "_R2_" | tr '\n' ' ')
    
    # Dry run: compact output
    if $dry_run; then
        echo "cat $R1_files > ${sample}_R1_concatenated.fastq.gz"
        echo "cat $R2_files > ${sample}_R2_concatenated.fastq.gz"
    else
        # Perform the actual concatenation
        echo "Concatenating files for sample $sample"
        cat $R1_files > ${sample}_R1_concatenated.fastq.gz
        cat $R2_files > ${sample}_R2_concatenated.fastq.gz
        echo "Done for sample $sample"
    fi
done

