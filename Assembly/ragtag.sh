#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=RAGTAG
#SBATCH --nodelist=node15
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=louis.dennu@obs-banyuls.fr
#SBATCH --mail-type=ALL

ASSEMBLY=$1
REF=$2
OUTPUT=$3

module load bioinfo/ragtag/2.1.0

ragtag.py scaffold -Cr -o $OUTPUT $REF $ASSEMBLY 

module purge

