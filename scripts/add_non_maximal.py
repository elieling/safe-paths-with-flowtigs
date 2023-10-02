
with open(snakemake.output.all_safe_paths, 'w') as destination:
    with open(snakemake.input.safe_paths, 'r') as source:
        for line in source:
            destination.write(line)
    if snakemake.wildcards.algorithm != "unitigs":
        with open(snakemake.input.unitigs, 'r') as source:
            for line in source:
                destination.write(line)




