#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=100000M
#SBATCH --nodes=28
#SBATCH --output=../to-tar/Human_gut

echo "Hello $USER! You are on node $HOSTNAME. The time is ${date}."

cp -r data/preprocessed_metagenome/Human_gut/ ../to-tar/
