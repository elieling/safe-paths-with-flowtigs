#!/usr/bin/env python3

from Bio import SeqIO
import numpy as np


# # contigs
# uni = 2168
# t_omni = 1235
# multi = 1240
# flow = 1461
# omni = 1219

# unaligned contigs
uni = 1314
t_omni = 135
multi = 137
flow = 241
omni = 133


unitigs = [i for i in SeqIO.parse("/home/ebingerv/VSC/snakemake-flowtigs/additional_scripts/data/unitigs.fasta", 'fasta')]
unitigs_arr = []
for sequence in unitigs:
    data = sequence.seq
    s_data = str(data)
    unitigs_arr.append(len(s_data))
unitigs_n = np.array(unitigs_arr)

gsum = -1
for i in range(0, 5000000):
    g = unitigs_n <= i
    gsum_early = gsum
    gsum = g.sum()
    if gsum == uni:
        print("exact", i, gsum)
    elif gsum > uni:
        break
print("unitigs", uni)
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
for i in range(0, 5000000):
    g = unitigs_n <= i
    gsum_early = gsum
    gsum = g.sum()
    if gsum == t_omni:
        print("exact", i, gsum)
    elif gsum > t_omni:
        break
print("trivial omnitigs", t_omni)
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
for i in range(0, 5000000):
    g = unitigs_n <= i
    gsum_early = gsum
    gsum = g.sum()
    if gsum == multi:
        print("exact", i, gsum)
    elif gsum > multi:
        break
print("muti-safe", multi)
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
for i in range(0, 5000000):
    g = unitigs_n <= i
    gsum_early = gsum
    gsum = g.sum()
    if gsum == flow:
        print("exact", i, gsum)
    elif gsum > flow:
        break
print("flowtigs", flow)
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
for i in range(0, 5000000):
    g = unitigs_n <= i
    gsum_early = gsum
    gsum = g.sum()
    if gsum == omni:
        print("exact", i, gsum)
    elif gsum > omni:
        break
print("omnitigs", omni)
print(i-1, gsum_early)
print(i, gsum)
print("#"*70)