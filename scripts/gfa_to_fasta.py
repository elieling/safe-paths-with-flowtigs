fasta_sequences = []
    
with open(snakemake.input.gfa, 'r') as gfa:
    for line in gfa:
        if line.startswith('S'):
            fields = line.strip().split('\t')
            sequence_id = fields[1]
            sequence = fields[2]
            fasta_sequences.append((sequence_id, sequence))

with open(snakemake.output.fasta, 'w') as fasta:
    for sequence_id, sequence in fasta_sequences:
        fasta.write(f'>{sequence_id}\n{sequence}\n')

