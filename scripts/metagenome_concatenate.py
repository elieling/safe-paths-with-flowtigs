#!/usr/bin/env python3

import pandas as pd

n_rows = snakemake.input.references.len()

for i in range(n_rows):
    with open(snakemake.input.references[i], 'r') as source:
        with open(snakemake.output.report, 'a') as destination:
            content = source.read()
            destination.write(content)