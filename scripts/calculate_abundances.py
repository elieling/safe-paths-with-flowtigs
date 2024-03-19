# Simulates non-zero abundances for a metagenome using the log-normal distribution

import os, shutil
import numpy as np
import math


genomes = []
folder = snakemake.input.reads
abundance_file = os.path.join(folder, "nanosim.abundances.tsv")
if os.path.exists(abundance_file):
    shutil.copy(abundance_file, snakemake.output.abundances)
    print("Abundance file found, so abundances are not simulated.")
else:
    for file in os.listdir(folder):
        if file.endswith(".fna") or file.endswith(".fasta") or file.endswith(".fq"):
            genomes.append(file)

    np.random.seed(0)
    results = []
    maximum = 0
    while len(results) < len(genomes):
        value = np.random.lognormal(0, 2)
        if value != 0.0: 
            results.append(math.ceil(value))
            maximum = max(maximum, math.ceil(value))
    print(results)
    genomes.sort()
    with open(snakemake.output.abundances, 'w') as outfile:
        outfile.write("Size\t50000\n")
        for i in range(len(genomes)):
            outfile.write(str(genomes[i]) + "\t" + str(results[i]) + "\n")

    print("Maximum:", maximum)