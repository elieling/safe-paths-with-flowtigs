import os
import numpy as np
import pandas as pd
import gzip
import shutil


counter = 1
for filename in os.listdir("/wrk-vakka/users/ebingerv/snakemake-flowtigs/data/meta/human_gut"):
    input_file = os.path.join("/wrk-vakka/users/ebingerv/snakemake-flowtigs/data/meta/human_gut/", filename)
    output_file = os.path.join("/wrk-vakka/users/ebingerv/snakemake-flowtigs/data/meta/Human_gut/", filename)
 
    if os.path.isdir(input_file): continue
    if os.path.isfile(output_file): continue

    with gzip.open(input_file, 'rb') as input:
        with open(output_file, 'wb') as output:
            shutil.copyfileobj(input, output)

    print(f"{counter}/898 done.")
    counter += 1


