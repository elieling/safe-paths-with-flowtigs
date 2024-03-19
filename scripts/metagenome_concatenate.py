#!/usr/bin/env python3

import pandas as pd
import os




# Collecting fasta files for metagenomes

abundances_df = pd.read_csv(snakemake.input.abundances, sep='\t')
names = abundances_df["Size"]
fasta_files = []
for name in names:
    print("#"*70)
    print(name)
    if name.endswith(".fna"):
        name = name[:-4] + ".fasta"
    if name.endswith(".fq"):
        name = name[:-3] + ".fasta"
    if name.endswith(".fasta"):
        fasta_files.append(os.path.join(snakemake.input.references, name))


print("¤"*70)
print(fasta_files)
# Concatenating the files

n_rows = len(fasta_files)

counter = 0
for i in range(n_rows):
    with open(fasta_files[i], 'r') as source:
        with open(snakemake.output.report, 'a') as destination:
            content = source.read()
            destination.write(content)
            print("Content: ", len(content))
            for character in content:
                if character in 'ACGT': counter += 1

with open(snakemake.output.number_of_characters, 'a') as number:
    number.write(str(counter))