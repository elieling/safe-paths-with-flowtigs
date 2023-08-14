#!/usr/bin/env python3

import pandas as pd
import os


# Collecting fasta files for metagenomes

abundance_df = pd.read_csv(snakemake.input.abundances, sep='\t')
names = abundance_df["Size"]
fasta_files = []
for name in names:
    if name.endswith(".fasta"):
        fasta_files.append(os.path.join(snakemake.input.references, name))



# Concatenating the files

n_rows = abundance_df.shape[0]

for i in range(n_rows):
    percentage = abundance_df.iloc[i, 1]
    if isinstance(percentage, int): n_occurences = percentage
    else: n_occurences = round(percentage)
    with open(fasta_files[i], 'r') as source:
        with open(snakemake.output.report, 'a') as destination:
            content = source.read()
            for j in range(n_occurences):
                destination.write(content)