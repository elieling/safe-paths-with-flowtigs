# Snakemake Flowtigs

A snakmake pipeline for comparing various safe paths with flowtigs in a De Bruijn graph of DNA reads in metagenomes. The pipeline also does all the pre-processing and the post-processing, and compares the results with structural contigs, extended contigs and unitigs.

## Input

The input dataset should be a folder located in the folder `data/meta/` from the root of this project. The input folder should contain all the genome reads as fasta files with the ending ".fasta" or ".fna". The folder can contain other files as long as they don't have these endings. The five input datasets used in the experiments can be found [here](https://zenodo.org/record/8434267).

## Output

The output of this pipeline is the LaTeX file `data/reports/output/meta_<name-of-the-input-folder>_k31ma1t28nm1/<report-name>/report.tex`. The figures of the latex file can be found in the folder `data/reports/hashdir`.

## Setup

The following has to be done once before running the pipeline.

- Install conda.
- Create the conda environment with the command `conda env create -f environment.yml` in the root of this project.

## Running Instructions

To run the pipeline:

- Activate the environment with the command `conda activate snakemake-flowtigs`.
- Execute the pipeline with the command `snakemake -j 1 "<path-to-this-project>/data/reports/output/meta_<name-of-the-input-folder>_k31ma1t28nm1/<report-name>/report.tex"`.

## Changing Parameters

The pipeline parameters can be changed by modifying the `k31ma1t28nm1` part of the running command. This will also change the name of the output file accordingly. The four parameters are represented by characters, followed by an integer indicating their value:
- k: size of the k-mers used in the De Bruijn graph.
- ma: minimum adundance.
- t: number of threads used.
- nm: binary value representing whether or not to append the unitigs to the results of the other safe-walk-algorithms. Use 1 to append unitigs, and 0 to not append them.
