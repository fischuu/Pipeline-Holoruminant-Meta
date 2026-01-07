#!/bin/bash
# Usage: ./split_bins.sh output_dir assembly1.fa.gz assembly2.fa.gz ...

set -euo pipefail  # exit on errors, undefined variables, or pipe failures
IFS=$'\n\t'

# First argument: output directory
out_dir="$1"
shift

# Remaining arguments: assemblies
assemblies=("$@")

mkdir -p "$out_dir"
echo "[INFO] Output directory: $out_dir"
echo "[INFO] Assemblies to process: ${assemblies[*]}"

# Loop over each assembly
for asm in "${assemblies[@]}"; do
    echo "[INFO] Processing assembly: $asm"

    # get assembly ID from the filename
    asm_id=$(basename "$asm" .fa.gz)
    echo "[INFO] Assembly ID: $asm_id"

    # Check if file exists
    if [[ ! -f "$asm" ]]; then
        echo "[ERROR] File not found: $asm"
        continue
    fi

    # Process the assembly
    gzip -dc "$asm" | awk -v outdir="$out_dir" -v asm="$asm_id" '
      BEGIN { 
        print "[INFO] Starting AWK processing for " asm > "/dev/stderr"
      }
      /^>/ {
        header=$0
        match(header, /bin_[0-9]+/)
        if (RSTART == 0) {
            print "[WARNING] No bin found in header: " header > "/dev/stderr"
            bin = "unknown"
        } else {
            bin = substr(header, RSTART, RLENGTH)
        }
        fname = asm ":" bin ".fa"
        print "[INFO] Writing header to file: " outdir "/" fname > "/dev/stderr"
        print header > outdir "/" fname
        next
      }
      {
        print $0 > outdir "/" fname
      }
      END {
        print "[INFO] Finished AWK processing for " asm > "/dev/stderr"
      }
    '
done

echo "[INFO] All assemblies processed."
