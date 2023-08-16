import os
import numpy as np
import pandas as pd
import requests



files_df = pd.read_csv("/wrk-vakka/users/ebingerv/snakemake-flowtigs/data/meta/human_gut_files.tsv", sep='\t')
accession_df = files_df["Assembly Accession"]
name_df = files_df["Assembly Name"]
rows = len(files_df.index)
for i in range(rows):
    accession = accession_df[i]
    name = name_df[i]
    if os.path.isfile(f"/wrk-vakka/users/ebingerv/snakemake-flowtigs/data/meta/human_gut/{accession}_{name}.fasta"): continue
    with open(f"/wrk-vakka/users/ebingerv/snakemake-flowtigs/data/meta/human_gut/{accession}_{name}.fasta", 'wb') as destination:
        url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{accession[4:7]}/{accession[7:10]}/{accession[10:13]}/{accession}_{name}/{accession}_{name}_genomic.fna.gz"
        destination.write(requests.get(url, allow_redirects=True).content)
    print(f"{i}/{rows} done: {url}")


