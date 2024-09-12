#!/usr/bin/env bash

# Determine the number of columns in the input file
num_columns=$(head -n1 unmapped.tsv | awk -F'\t' '{print NF}')

# Loop through each column and extract the corresponding elements
for ((col=1; col <= num_columns; col++)); do
# Extract the elements from the current column, split by '_' and join with tab delimiter
awk -F'\t' -v col="$col" '
    NR==1 {
    split($col, arr, "_");
    printf("%s\t%s\t%s\n", arr[1], arr[2], arr[3]);
    }
    NR>1 {
    value = 100 - $col;
    printf("%.5f\n", value);
    }
' unmapped.tsv | paste -d'\t' -s
done > longer.tsv
