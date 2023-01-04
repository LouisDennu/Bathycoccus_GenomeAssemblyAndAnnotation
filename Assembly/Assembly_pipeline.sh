#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=Bathycoccus_genome_assembly
#SBATCH --nodelist=node24
#SBATCH --cpus-per-task=8
#SBATCH --mail-user='dennu.louis@gmail.com'
#SBATCH --mail-type=ALL

NAME=$1
FLYE=$2
PILON=$3

NANOPORE_READS=$4
ILLUMINA_READS_1=$5
ILLUMINA_READS_2=$6
REFERENCE=$7

mkdir /scratch/ldennu
mkdir /scratch/ldennu/Bathycoccus_assembly
mkdir /scratch/ldennu/Bathycoccus_assembly/${NAME}

cd /scratch/ldennu/Bathycoccus_assembly/${NAME}

##########
## FLYE assembly of nanopore reads
time python ${FLYE} --nano-hq ${NANOPORE_READS} --genome-size 15m --scaffold --threads 8 --out-dir FLYE

##########
## MEDAKA correction with nanopore reads
module load bioinfo/medaka/1.5

time medaka_consensus -i ${NANOPORE_READS} -d ./FLYE/assembly.fasta -o MEDAKA m r941_min_high_g360 -t 2

##########
## PILON polishing with illumina reads
module load bioinfo/bwa/0.7.17
module load bioinfo/samtools/1.9

bwa index ./MEDAKA/consensus.fasta
bwa mem   ./MEDAKA/consensus.fasta ${ILLUMINA_READS_1} ${ILLUMINA_READS_2} | samtools view -b | samtools sort -o ${NAME}_FLYE_MEDAKA.sorted.bam
samtools index ${NAME}_FLYE_MEDAKA.sorted.bam
java -Xmx16G -jar ${PILON} --genome ./MEDAKA/consensus.fasta --frags ${NAME}_FLYE_MEDAKA.sorted.bam --outdir PILON

##########
##RAGTAG scaffolding with reference assembly
module load bioinfo/ragtag/2.1.0

ragtag.py scaffold -Cr -o RAGTAG ${REFERENCE} ./PILON/pilon.fasta

##########
module purge
