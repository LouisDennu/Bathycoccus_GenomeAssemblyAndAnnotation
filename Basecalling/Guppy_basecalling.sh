#!/bin/bash
#SBATCH -J bonito 
#SBATCH -p gpu
#SBATCH -A gpu_group
#SBATCH -c 8
INPUT=$1
OUTPUT=$2
CUDA=$3

#loading modules
module load bioinfo/guppy-gpu/6.1.2
module load bioinfo/nanoplot/1.19.0

#running basecalling
time guppy_basecaller -c dna_r9.4.1_450bps_sup.cfg -i ${INPUT} -r -s ${OUTPUT} --barcode_kits "SQK-RBK004" --compress_fastq --num_callers 8 --gpu_runners_per_device 8 --min_qscore 7 -x cuda:${CUDA}
cd ${OUTPUT}
NanoPlot -t 4 -o ./nanoplot --summary ./sequencing_summary.txt
