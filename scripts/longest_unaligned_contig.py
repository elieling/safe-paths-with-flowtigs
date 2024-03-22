#!/usr/bin/env python3

import pandas as pd
import sys, os


report = sys.argv[1]
unaligned = sys.argv[2]
output_file_name = sys.argv[3]

if os.path.getsize(unaligned) > 0:
    df = pd.read_csv(unaligned, delimiter='\t', engine='python')
    print(df)
    unaligned_length = df["Unaligned_length"]
    greatest_value = unaligned_length.max()
else: greatest_value = 0
with open(report, 'r') as source:
    with open(output_file_name, 'a') as destination:
        for line in source:
            if line.startswith("Genome fraction"):
                destination.write(f"Longest unaligned contig & {greatest_value} \\\\ \hline\n")
            destination.write(line)

#output_file = open(output_file_name, 'w')