#!/bin/bash -l
#SBATCH --time=00:10:00
#SBATCH --mem=10G

echo "Hello $USER! You are on node $HOSTNAME. The time is ${date}. Ready to extract files!"

srun python3 scripts/extract_human_gut_files.py
