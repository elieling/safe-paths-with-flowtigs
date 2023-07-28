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

gsum = -1
for i in range(30, 5000000):
    g = unitigs_n > i
    gsum_early = gsum
    gsum = g.sum()
    if gsum == 1314:
        print(i)
    elif gsum < 1314:
        break
print("unitigs 1314")
print(i-1, gsum_early)
print(i, gsum)
print("#"*70)


unitigs = [i for i in SeqIO.parse("/home/ebingerv/VSC/snakemake-flowtigs/additional_scripts/data/flowtigs.fasta", 'fasta')]
unitigs_arr = []
for sequence in unitigs:
    data = sequence.seq
    s_data = str(data)
    unitigs_arr.append(len(s_data))
unitigs_n = np.array(unitigs_arr)

gsum = -1
for i in range(30, 5000000):
    g = unitigs_n > i
    gsum_early = gsum
    gsum = g.sum()
    if gsum == 241:
        print("exact", i, gsum)
    elif gsum < 241:
        break
print("flowtigs 241")
print(i-1, gsum_early)
print(i, gsum)
print("#"*70)


unitigs = [i for i in SeqIO.parse("/home/ebingerv/VSC/snakemake-flowtigs/additional_scripts/data/trivial-omnitigs.fasta", 'fasta')]
unitigs_arr = []
for sequence in unitigs:
    data = sequence.seq
    s_data = str(data)
    unitigs_arr.append(len(s_data))
unitigs_n = np.array(unitigs_arr)

gsum = -1
for i in range(30, 5000000):
    g = unitigs_n > i
    gsum_early = gsum
    gsum = g.sum()
    if gsum == 135:
        print("exact", i, gsum)
    elif gsum < 135:
        break
print("trivial omnitigs 135")
print(i-1, gsum_early)
print(i, gsum)
print("#"*70)


unitigs = [i for i in SeqIO.parse("/home/ebingerv/VSC/snakemake-flowtigs/additional_scripts/data/multi-safe.fasta", 'fasta')]
unitigs_arr = []
for sequence in unitigs:
    data = sequence.seq
    s_data = str(data)
    unitigs_arr.append(len(s_data))
unitigs_n = np.array(unitigs_arr)

gsum = -1
for i in range(30, 5000000):
    g = unitigs_n > i
    gsum_early = gsum
    gsum = g.sum()
    if gsum == 137:
        print("exact", i, gsum)
    elif gsum < 137:
        break
print("muti-safe 137")
print(i-1, gsum_early)
print(i, gsum)
print("#"*70)


unitigs = [i for i in SeqIO.parse("/home/ebingerv/VSC/snakemake-flowtigs/additional_scripts/data/omnitigs.fasta", 'fasta')]
unitigs_arr = []
for sequence in unitigs:
    data = sequence.seq
    s_data = str(data)
    unitigs_arr.append(len(s_data))
unitigs_n = np.array(unitigs_arr)

gsum = -1
for i in range(30, 5000000):
    g = unitigs_n > i
    gsum_early = gsum
    gsum = g.sum()
    if gsum == 133:
        print("exact", i, gsum)
    elif gsum < 133:
        break
print("omnitigs 133")
print(i-1, gsum_early)
print(i, gsum)
print("#"*70)