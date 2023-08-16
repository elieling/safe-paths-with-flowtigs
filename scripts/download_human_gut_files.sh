#!/bin/bash -l
#SBATCH --time=00:10:00
#SBATCH --mem=100M

echo "Hello $USER! You are on node $HOSTNAME. The time is ${date}. Ready to download files!"

srun python3 scripts/download_human_gut_files.py
