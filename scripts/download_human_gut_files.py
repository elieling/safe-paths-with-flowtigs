import os
import numpy as np
import pandas as pd
import requests
from joblib import Parallel, delayed


def download_single_file(accession, name):
    if not os.path.isfile(os.path.join(snakemake.params.path, f"/human_gut/{accession}_{name}.fasta")): 
        with open(os.path.join(snakemake.params.path, f"human_gut/{accession}_{name}.fasta"), 'wb') as destination:
            url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{accession[4:7]}/{accession[7:10]}/{accession[10:13]}/{accession}_{name}/{accession}_{name}_genomic.fna.gz"
            destination.write(requests.get(url, allow_redirects=True).content)


os.mkdir(os.path.join(snakemake.params.path, "human_gut"))
files_df = pd.read_csv(os.path.join(snakemake.params.path, "human_gut_files.tsv"), sep='\t')
accession_df = files_df["Assembly Accession"]
name_df = files_df["Assembly Name"]
rows = len(files_df.index)
Parallel(n_jobs=28)(delayed(download_single_file)(accession_df[i], name_df[i]) for i in range(rows))


