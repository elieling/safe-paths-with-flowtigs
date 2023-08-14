# Simulates non-zero abundances for a metagenome using the Poisson distribution

import os
import numpy as np


genomes = []
folder = snakemake.input.reads
for file in os.listdir(folder):
    if file.endswith(".fna") or file.endswith(".fasta"):
        genomes.append(file)

np.random.seed(42)
# poisson_distribution = np.random.poisson(1, len(genomes))
# print(poisson_distribution)
results = []
while len(results) < len(genomes):
    value = np.random.poisson(1)
    if value != 0: results.append(value)
print(results)
with open(snakemake.output.abundances, 'w') as outfile:
    # outfile.write("Size\t50000\n")
    # for i in range(len(genomes)):
    #     outfile.write(str(genomes[i]) + "\t" + str(poisson_distribution[i]) + "\n")
    outfile.write("Size\t50000\n")
    for i in range(len(genomes)):
        outfile.write(str(genomes[i]) + "\t" + str(results[i]) + "\n")