#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np


unitigs = [i for i in SeqIO.parse("/home/ebingerv/VSC/snakemake-flowtigs/additional_scripts/data/unitigs.fasta", 'fasta')]
unitigs_arr = []
for sequence in unitigs:
    data = sequence.seq
    s_data = str(data)
    unitigs_arr.append(len(s_data))
unitigs_n = np.array(unitigs_arr)

for i in range(30, 5000000):
    print(i)
    g = unitigs_n > i
    gsum = g.sum()
    print(gsum)
    if gsum == 1314:
        print(i)
    elif gsum < 1314:
        break
    if (1==100):
        print(i, gsum)
    if (1==1000):
        print(i, gsum)
    if (1==10000):
        print(i, gsum)
