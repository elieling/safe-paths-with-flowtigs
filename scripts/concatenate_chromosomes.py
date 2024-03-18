#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import os

print("¤"*100)
print(snakemake.output.report)
os.makedirs(snakemake.output.report)
if os.path.exists(snakemake.output.report): print("YES")
else: print("NO")

files = os.listdir(snakemake.input.assembly)
fasta_files = [file for file in files if file.endswith(".fna") or file.endswith(".fasta")]

for file in fasta_files:
    filename = os.path.basename(file)
    with open(os.path.join(snakemake.input.assembly, filename), 'r') as infile:
        with open(os.path.join(snakemake.output.report, filename), 'w') as outfile:

            outfile.write('>Sequence\n')

            # Going through all the sequences from the input file
            # sequences=[i for i in SeqIO.parse(infile, 'fasta')]
            for sequence in SeqIO.parse(infile, 'fasta'):

                # Extracting the actual sequence and transforming it to string format
                data = sequence.seq
                s_data = str(data)

                # Writing our result to the output file
                outfile.write('{}'.format(Seq(s_data)))

            outfile.write('\n')
        
 



