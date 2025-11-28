#!/bin/bash
#SBATCH --job-name=bathycoccus_assembly_master
#SBATCH --partition=normal
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/cluster/master.out
#SBATCH --error=logs/cluster/master.err

# Load necessary modules
module load singularity/4.0.1
module load python/3.8.12

# Launch Snakemake with the Slurm profile
snakemake --profile profiles/slurm -j 24
