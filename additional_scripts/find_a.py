#!/usr/bin/env python3

from Bio import SeqIO
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
    if gsum == 2168:
        print(i)
    elif gsum < 2168:
        break
print("unitigs 2168")
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
    if gsum == 1461:
        print("exact", i, gsum)
    elif gsum < 1461:
        break
print("flowtigs 1461")
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
    if gsum == 1235:
        print("exact", i, gsum)
    elif gsum < 1235:
        break
print("trivial omnitigs 1235")
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
    if gsum == 1240:
        print("exact", i, gsum)
    elif gsum < 1240:
        break
print("muti-safe 1240")
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
    if gsum == 1219:
        print("exact", i, gsum)
    elif gsum < 1219:
        break
print("omnitigs 1219")
print(i-1, gsum_early)
print(i, gsum)
print("#"*70)