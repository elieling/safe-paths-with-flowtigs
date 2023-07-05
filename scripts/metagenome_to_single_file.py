#!/usr/bin/env python3

import pandas as pd


abundance_df = pd.read_csv(snakemake.input.abundances, sep='\t')
n_rows = abundance_df.shape[0]

for i in range(n_rows):
    percentage = abundance_df.iloc[i]
    n_occurences = round(percentage)
    with open(snakemake.input.references[i], 'r') as source:
        with open(snakemake.output.report, 'a') as destination:
            content = source.read()
            for j in range(n_occurences):
                destination.write(content)

