# Experiment pipeline

This is an experiment pipeline to compare the performance of [flowtigs](https://www.biorxiv.org/content/10.1101/2023.11.17.567499v1) with other safety algorithms.

## Setting up the environment

To set up the environment, download [conda](https://www.anaconda.com/), then in the project directory, create the environment with
```
conda env create environment.yml
```
then activate the environment with
```
conda activate snakemake-flowtigs
```

## Reproducing experiments with simulated data

First, ensure you are using the correct version by running in the project directory
```
git checkout a2698923762b58129fd027602f242209cecdd2d7
```
Then, include the reference of the desired metagenome in the folder
```
safe-paths-with-flowtigs/data/meta/<name_of_the_dataset>
```
The folder should contain a file ending in `.fna` or `.fasta` for each genome in the reference metagenome, as well as a file named `nanosim.abundances.tsv` containing the abundance information in tsv format, so that the first column contains the name of the genome and the second column contains its abundance in percentage. If the abundance file is missing, abundances will be simulated. The datasets used for the experiments can be found [here](https://zenodo.org/records/8434267). Then, the file 
```
safe-paths-with-flowtigs/data/meta/<name_of_the_dataset>
```
should be copied in 
```
safe-paths-with-flowtigs/data/meta/<name_of_the_dataset>_reference
```
because we also simulate the reads when working with simulated data. Then, the pipeline is run with the command

```
snakemake --use-conda -j 1 "<path>/safe-paths-with-flowtigs/data/reports/output/meta_<metagenome>_k<k>ma1t<threads>nm1/<name_of_the_report>/report.tex"
```

The output report will be located in the file 
```
data/reports/output/meta_<metagenome>\_k<k>ma1t<threads>nm1/<name_of_the_report>/report.tex
```

## Reproducing experiments with real data

To ensure you are using the correct version, run in the project directory
```
git checkout bd04d0789f313415293ac3bdace399c6b0c035e5
```
Then, include the desired reads in the file
```
safe-paths-with-flowtigs/data/meta/<name_of_the_dataset>/reads.fq
```
 and the reference of the desired metagenome in the folder
 ```
safe-paths-with-flowtigs/data/meta/<name_of_the_dataset>_reference
```
The reference metagenome should be in the same format as with simulated data, but does not need an abundance file. In our experiments, we use the [Zymo D6331](https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards/products/zymobiomics-gut-microbiome-standard) standard gut microbiome and reads with accession [SRR13128014](https://www.ebi.ac.uk/ena/browser/view/SRX9569057) and we use a minimum abundance threshold of 5, i.e., we exclude all edges that have a multiplicity lower than 5. Then, the pipeline is run with the command

```
snakemake --use-conda -j 1 "<path>/safe-paths-with-flowtigs/data/reports/output/real_<metagenome>_k<k>ma<minimum_abundance>t<threads>nm0th<threshold>/<name_of_the_report>/report.tex" 
```

The output report will be located in the file 
```
data/reports/output/meta_<metagenome>_k<k>ma<minimum_abundance>t<threads>nm0th<threshold>/<name_of_the_report>/report.tex
```

To run the pipeline without using the additionnal filtering, instead run

```
snakemake --use-conda -j 1 "<path>/safe-paths-with-flowtigs/data/reports/output_no_filtering/real_<metagenome>_k<k>ma<minimum_abundance>t<threads>nm0th<threshold>/<name_of_the_report>/report.tex" 
```

The output report will then be located in the file 
```
data/reports/output_no_filtering/meta_<metagenome>_k<k>ma<minimum_abundance>t<threads>nm0th<threshold>/<name_of_the_report>/report.tex
```
