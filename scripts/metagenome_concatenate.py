#!/usr/bin/env python3

import pandas as pd
import os




# Collecting fasta files for metagenomes

abundances_df = pd.read_csv(snakemake.input.abundances, sep='\t')
names = abundances_df["Size"]
fasta_files = []
for name in names:
    if name.endswith(".fasta"):
        fasta_files.append(os.path.join(snakemake.input.references, name))



# Concatenating the files

n_rows = len(fasta_files)

for i in range(n_rows):
    with open(fasta_files[i], 'r') as source:
        with open(snakemake.output.report, 'a') as destination:
            content = source.read()
            destination.write(content)