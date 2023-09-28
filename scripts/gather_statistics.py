import pandas as pd
import re

statistics = []
with open(snakemake.input.log, 'r') as input_file:
        for line in input_file:                
            if re.match(r'.*Graph has \d+ nodes and \d+ edges', line):
                line = line[line.find("Graph"):]
                statistics = re.findall(r'\d+', line)
                break

        assert "statistics" in locals(), f"No statistics found in {snakemake.input.log}"



with open(snakemake.input.edges_in_cycles, 'r') as input_file:
        for line in input_file:                
            if re.match(r'.*Cycles contain a total of \d+ edges', line):
                line = line[line.find("Cycles"):]
                statistics.append(re.search(r'\d+', line).group(0))
                break



df = pd.DataFrame({'nodes': [statistics[0]], 'edges': [statistics[1]], 'edges in cycles': [statistics[2]]})
df.to_csv(snakemake.output.statistics, sep='\t')