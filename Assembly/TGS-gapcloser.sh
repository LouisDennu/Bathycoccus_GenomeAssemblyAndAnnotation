#!/bin/bash
#SBATCH -p normal
#SBATCH --job-name=TGS-gapcloser
#SBATCH --nodelist=node15
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=louis.dennu@obs-banyuls.fr
#SBATCH --mail-type=ALL

ASSEMBLY=$1
READS=$2
OUTPUT=$3
ILLMUNINA1=$4
ILLUMINA2=$5

#TGS-gapcloser_path="/scratch/ldennu/gapcloser/TGS-GapCloser/TGS-GapCloser.sh" 
#pilon_path="/scratch/ldennu/pilon/pilon-1.24.jar"
#java_path="/usr/java/default/bin/java"
#samtools_path="/usr/local/samtools-1.9/bin/samtools"

module load bioinfo/seqtk/1.3-r106

mkdir /scratch/ldennu/gapcloser/$(basename ${ASSEMBLY})_GC
cd /scratch/ldennu/gapcloser/$(basename ${ASSEMBLY})_GC
rsync -vaurL $READS .

seqtk seq -a $READS > ONTreads.fasta
${TGS-gapcloser_path} --scaff $ASSEMBLY --reads ./ONTreads.fasta --output $3 --pilon $pilon_path --javaj $java_path --ngs $4 $5 --samtools $samtools_path --p_round 1 --r_round 0 --min_nread 3 --pilon_meme 16G

module purge

