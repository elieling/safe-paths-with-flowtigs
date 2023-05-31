#!/usr/bin/env python3

from Bio import SeqIO
# from StringIO import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

        
k_string = snakemake.wildcards.k
k = int(k_string)
record_list = []
prints = []
print("Hello World")
prints.append("Hello World in terminal")
with open(snakemake.input.assembly, 'r') as infile:
    sequences=[i for i in SeqIO.parse(infile, 'fasta')]
    s_record = str(sequences[0])
    first_elements = s_record[:(k-1)]
    s_record = s_record + first_elements
    print("#"*70,"\ns_recors:",s_record,"\nlen(sequences):",len(sequences),"\nsequences[0]:",sequences[0].name, "\n","#"*70)
    print("Nmae:",sequences[0].name)
    data = sequences[0].seq
    print("### Sequence:", data)
    s_data = str(data)
    print("### Sequence_s:", s_data)
    s_data = s_data + 'A'
    print("### Sequence_s:", s_data)
    sequences[0].seq = s_data
    data = sequences[0].seq
    print("### Sequence_s:", data)
    # print("### Sequence_s:", sequences[0])
    for record in SeqIO.parse(infile, 'fasta'):
        # transform "record" from Bio.SeqRecord.SeqRecord to string
        s_record = str(record)
        if (len(s_record) < k - 1):
            print("Error: Sequence is too short compared to k.")
            prints.append("Error: Sequence is too short compared to k.")
            record_list.append(Seq("Error: Sequence is too short compared to k."))
            # record_list.append(SeqIO.read("Error: Sequence is too short compared to k.", "fasta"))
        else:
            # Copying the first k-1 elements of the record to the end to make it circular
            first_elements = s_record[:(k-1)]
            s_record = s_record + first_elements
            print("#"*70,"\n",s_record,len(sequences),sequences[0].name, "\n","#"*70)
            prints.append(s_record)
            prints.append(len(sequences))
            prints.append(sequences[0].name)
            record_list.append(Seq(s_record))


with open(snakemake.output.report, 'w') as outfile:
    outfile.write('contig_count: {}\n'.format(data))
    # SeqIO.write(sequences[0],outfile)
    # outfile.write('contig_count: {}\n'.format(record_list[0]))


    # Kind of prints the sequence in the fasta file, find out how to extract the interesting info. Google.
    # When done, write documentation, run entire pipeline, see results and report to Zeb whether or not circular. Then push to github.