import os
import numpy as np
import pandas as pd
import random
import shutil


n_files = int(snakemake.wildcards.n_files)

os.mkdir(snakemake.output.selected_files)
files_df = pd.read_csv(os.path.join(snakemake.input.file_list), sep='\t')
accession_df = files_df["Assembly Accession"]
name_df = files_df["Assembly Name"]

random.seed(0)
random_file_indexes = [random.randint(0, 897) for _ in range(n_files)]

for file_index in random_file_indexes:
    accession = accession_df[file_index]
    name = name_df[file_index]
    file_name = f"{accession}_{name}.fasta"
    shutil.copy(os.path.join(snakemake.input.all_files, file_name), snakemake.output.selected_files)