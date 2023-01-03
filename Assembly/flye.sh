#!/bin/bash
#SBATCH --partition=normal      ### Partition
#SBATCH --job-name=FLYE    ### Job Name
#SBATCH --nodelist=node15
#SBATCH --cpus-per-task=12
#SBATCH --mail-user='dennu.louis@gmail.com'
#SBATCH --mail-type=ALL

READS=$1
OUTDIR=$2

for reads in $(ls ${READS}* | grep -e ".fastq.gz$")
do
	name=$(basename $reads .fastq.gz)
	mkdir ${OUTDIR}${name}_FLYE
	python /scratch/ldennu/flye/Flye/bin/flye --nano-hq $reads --genome-size 15m --scaffold --threads 12 --out-dir ${OUTDIR}${name}_FLYE
	echo "$name done"
done
module purge

