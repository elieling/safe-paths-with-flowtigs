#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq

# Size of k in int format     
k_string = snakemake.wildcards.k
k = int(k_string)

record_list = []
with open(snakemake.input.assembly, 'r') as infile:

    # Creating a list of all the sequences from the input file
    sequences=[i for i in SeqIO.parse(infile, 'fasta')]
    for sequence in sequences:

        # Extracting the actual sequence and transforming it to string format
        data = sequence.seq
        s_data = str(data)

        # The circularization works only if the size of the sequence s at least is at least k-1.
        if (len(s_data) < k - 1):
            print("Error: Sequence is too short compared to k.")
            record_list.append(Seq("Error: Sequence is too short compared to k."))
        else:
            # Copying the first k-1 elements of the record to its end to make it circular
            first_elements = s_data[:(k-1)]
            s_data = s_data + first_elements

            # sequence.seq = s_data
            record_list.append(Seq(s_data))


# Writing our result to the output file
with open(snakemake.output.report, 'w') as outfile:
    for i in range(len(record_list)):
        outfile.write('>Sequence{}\n'.format(i))
        outfile.write('{}\n'.format(record_list[i]))
 