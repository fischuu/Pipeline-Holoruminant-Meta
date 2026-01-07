#!/bin/bash
# Usage: ./split_bins.sh output_dir assembly1.fa.gz assembly2.fa.gz ...

# First argument: output directory
out_dir="$1"
shift

# Remaining arguments: assemblies
assemblies=("$@")

mkdir -p "$out_dir"

# Loop over each assembly
for asm in "${assemblies[@]}"; do
    # get assembly ID from the filename
    asm_id=$(basename "$asm" .fa.gz)

    gzip -dc "$asm" | awk -v outdir="$out_dir" -v asm="$asm_id" '
      /^>/ {
        header=$0
        # extract bin_XXXXX from header
        match(header, /bin_[0-9]+/, b)
        bin = b[0]
        fname = asm ":" bin ".fa"
        print header > outdir "/" fname
        next
      }
      {
        print $0 > outdir "/" fname
      }
    '
done
