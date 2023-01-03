#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=MEDAKA
#SBATCH --nodelist=node15
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=louis.dennu@obs-banyuls.fr
#SBATCH --mail-type=ALL

ASSEMBLY=$1
FASTQ=$2
OUTDIR=$3

module load bioinfo/medaka/1.5

time medaka_consensus -i ${FASTQ} -d ${ASSEMBLY} -o ${OUTDIR} -m r941_min_high_g360 -t 2

module purge

