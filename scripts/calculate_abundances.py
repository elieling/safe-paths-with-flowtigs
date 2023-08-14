# Simulates non-zero abundances for a metagenome using the log-normal distribution

import os
import numpy as np


genomes = []
folder = snakemake.input.reads
for file in os.listdir(folder):
    if file.endswith(".fna") or file.endswith(".fasta"):
        genomes.append(file)

np.random.seed(42)
results = []
while len(results) < len(genomes):
    value = np.random.lognormal(0, 1)
    if value != 0.0: results.append(value)
print(results)
with open(snakemake.output.abundances, 'w') as outfile:
    outfile.write("Size\t50000\n")
    for i in range(len(genomes)):
        outfile.write(str(genomes[i]) + "\t" + str(results[i]) + "\n")