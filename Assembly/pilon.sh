#!/bin/bash
#SBATCH --partition=normal      ### Partition
#SBATCH --job-name=PILON    ### Job Name
#SBATCH --nodelist=node15
#SBATCH --cpus-per-task=4
#SBATCH --mail-user='dennu.louis@gmail.com'
#SBATCH --mail-type=ALL

module load bioinfo/bwa/0.7.17
module load bioinfo/samtools/1.9

ASSEMBLY=$1
READS=$2
OUTDIR=$3

for assembly in $(ls ${ASSEMBLY}* | grep -e ".fasta$")
do
	name=$(basename $assembly .fasta)
	mkdir ${OUTDIR}${name}_PILON
	cd ${OUTDIR}${name}_PILON
	bwa index $assembly
	bwa mem $assembly ${READS}$(echo $name | grep "^...." -o -m 1)* | samtools view -b | samtools sort -o $name.sorted.bam
	samtools index $name.sorted.bam
	java -Xmx16G -jar /scratch/ldennu/pilon/pilon-1.24.jar --genome $assembly --frags $name.sorted.bam --outdir ${OUTDIR}${name}_PILON
	echo "$name done"
done
module purge

