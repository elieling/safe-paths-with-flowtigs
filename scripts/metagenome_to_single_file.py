#!/usr/bin/env python3

import pandas as pd


record_list = []
abundance_df = pd.read_csv(snakemake.input.abundances, sep='\t')
n_rows = abundance_df.shape[0]
print(n_rows)

for i in range(n_rows):
    print("I", i)
    print(abundance_df)
    print(abundance_df.iloc[0,0])
    print(abundance_df.iloc[i])
    percentage = abundance_df.iloc[i, 0]
    n_occurences = round(percentage)
    with open(snakemake.input.references[i], 'r') as source:
        with open(snakemake.output.report, 'a') as destination:
            content = source.read()
            for j in range(n_occurences):
                destination.write(content)


# with open(snakemake.input.abundances, 'r') as infile:

#     # Going through all the sequences from the input file
#     sequences=[i for i in SeqIO.parse(infile, 'fasta')]
#     for sequence in sequences:

#         # Extracting the actual sequence and transforming it to string format
#         data = sequence.seq
#         s_data = str(data)

#         # The circularization works only if the size of the sequence s at least is at least k-1.
#         if (len(s_data) < k - 1):
#             print("Error: Sequence is too short compared to k.")
#             record_list.append(Seq("Error: Sequence is too short compared to k."))
#         else:
#             # Copying the first k-1 elements of the record to its end to make it circular
#             first_elements = s_data[:(k-1)]
#             s_data = s_data + first_elements

#             record_list.append(Seq(s_data))


# # Writing our result to the output file
# with open(snakemake.output.report, 'w') as outfile:
#     for i in range(len(record_list)):
#         outfile.write('>Sequence{}\n'.format(i))
#         outfile.write('{}\n'.format(record_list[i]))
 