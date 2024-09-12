#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
import sys

input_fasta = sys.argv[1]
input_mapping = sys.argv[2]

fasta = {x.id: str(x.seq) for x in SeqIO.parse(input_fasta, format="fasta")}

name_mapping = pd.read_table(input_mapping)
name_mapping_dict = {
    seq_id: sequence
    for seq_id, sequence in name_mapping[["contig", "seqname"]].values.tolist()
}

new_fasta = (
    f">{name_mapping_dict[old_identifier]}\n{sequence}"
    for old_identifier, sequence in fasta.items()
    if old_identifier in name_mapping_dict.keys()
)

for record in new_fasta:
    sys.stdout.write(record + "\n")
