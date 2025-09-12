#!/usr/bin/env python3
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as fin, open(output_file, "w") as fout:
    header = fin.readline()
    fout.write(header)
    for line in fin:
        parts = line.rstrip("\n").split("\t")
        if len(parts) > 1:
            parts[1] = parts[0].split("@")[0]  # scaffold = alles vor "@"
        fout.write("\t".join(parts) + "\n")
