#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq

# Size of k in int format     
k_string = snakemake.wildcards.k
k = int(k_string)
counter = 0

record_list = []
with open(snakemake.input.assembly, 'r') as infile:
    with open(snakemake.output.report, 'w') as outfile:

        # Going through all the sequences from the input file
        # sequences=[i for i in SeqIO.parse(infile, 'fasta')]
        for sequence in SeqIO.parse(infile, 'fasta'):

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

                # record_list.append(Seq(s_data))
            # Writing our result to the output file
            outfile.write('>Sequence{}\n'.format(counter))
            outfile.write('{}\n'.format(Seq(s_data)))
            counter += 1

        
 
