import os
import numpy as np
import pandas as pd
import gzip
import shutil
from joblib import Parallel, delayed


def extract_single_file(filename):
    input_file = os.path.join(snakemake.params.path, "human_gut/", filename)
    output_file = os.path.join(snakemake.params.path, "Human_gut/", filename)
 
    if os.path.isdir(input_file): return
    if os.path.isfile(output_file): return

    with gzip.open(input_file, 'rb') as input:
        with open(output_file, 'wb') as output:
            shutil.copyfileobj(input, output)


os.mkdir(os.path.join(snakemake.params.path, "Human_gut"))
Parallel(n_jobs=28)(delayed(extract_single_file)(filename) for filename in os.listdir(os.path.join(snakemake.params.path, "human_gut")))

     

