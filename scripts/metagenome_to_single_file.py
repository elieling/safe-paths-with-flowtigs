#!/usr/bin/env python3

import pandas as pd
import logging


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


abundance_df = pd.read_csv(snakemake.input.abundances, sep='\t')
logger.info(abundance_df)
logger.info("len %s", len(abundance_df))
n_rows = abundance_df.shape[0]

for i in range(n_rows):
    percentage = abundance_df.iloc[i]
    n_occurences = round(percentage)
    with open(snakemake.input.references[i], 'r') as source:
        with open(snakemake.output.report, 'a') as destination:
            content = source.read()
            for j in range(n_occurences):
                destination.write(content)

