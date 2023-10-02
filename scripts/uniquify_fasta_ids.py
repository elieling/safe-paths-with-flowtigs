import sys
from Bio import SeqIO



input_file = snakemake.input.safe_paths
output_file = snakemake.output.uniquified

def read_fasta(input_file):
    id = 0
    for record in SeqIO.parse(input_file, "fasta"):
        record.id = str(id) + "_" + record.id
        id += 1
        yield record

SeqIO.write(read_fasta(input_file), output_file, "fasta")