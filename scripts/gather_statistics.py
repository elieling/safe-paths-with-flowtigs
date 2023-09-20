import pandas as pd

with open(snakemake.input, 'r') as input_file:
        for line in input_file:
            if re.match(r'Graph has \d+ nodes and \d+ edges', line):
                line = line[line.find("Graph"):]
                statistics = re.finall(r'\d+', line)
                break

        assert "statistics" in locals(), f"No statistics found in {snakemake.input}"


df = pd.DataFrame({'nodes': [statistics[0]], 'edges': [statistics[1]]})
df.to_csv(snakemake.output, sep='\t')