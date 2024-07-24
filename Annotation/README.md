# Bathycoccus genome annotation

## Model training and structural annotation

<p align="center">
<img src="https://github.com/LouisDennu/BathycoccusPangenome/blob/main/Annotation/Bathycoccus_pangenome_pipelines-Genome_annotation_pipeline.png">
</p>

### Repeat masking

Masking of repeat sequences was performed using RepeatMasker 4.1.2 (https://github.com/rmhubley/RepeatMasker) and a non-redundant repeat library created from individual genomes.

```bash
RepeatMasker -xsmall -lib {input.repeat_db} -s -pa 12 {input.assembly} --dir /scratch/ldennu/assembly
```

### Transcript assembly

Available Illumina RNAseq whole transcriptome data of *B. prasinos* strain RCC4752 (https://www.ncbi.nlm.nih.gov/sra/SRX554258[accn]) were mapped upon RCC4752 genome assembly using HISAT2 2.1.0 (http://daehwankimlab.github.io/hisat2/) and Samtools 1.9 (https://github.com/samtools/samtools).

```bash
hisat2-build -p 12 {input.soft_masked_assembly} /scratch/ldennu/RNAseq_snakemake/{wildcards.strain}_masked.index

hisat2 -x /scratch/ldennu/RNAseq_snakemake/{wildcards.strain}_masked.index -1 {input.RNAseq_illumina_1} -2 {input.RNAseq_illumina_2} -S {output.sam_file} -p 12

samtools view -@ 12 -Sb {output.sam_file} | samtools sort -o {output.bam_file}
```

Transcript assembly was performed using Stringtie 1.3.4 (https://github.com/gpertea/stringtie)

```bash
stringtie {input.bam_file} -p 12 -v -o {output.transcripts_gft}
```

### Augustus and Genemark-ETP model training using BRAKER3

GENEMARK-ETP 4.71 (https://github.com/gatech-genemark/GeneMark-ETP) and AUGUSTUS 3.3.3 (https://github.com/Gaius-Augustus/Augustus) gene prediction models were trained through a run of BRAKER3 2.1.6 (https://github.com/Gaius-Augustus/BRAKER) using assembled transcripts and plant protein database (plant_odb10) as evidences.

```bash
braker.pl \
		--genome {input.assembly} \
		--bam {input.bam} \
		--prot_seq {input.protein_db} \
		--etpmode \
		--cores 12 \
		--GENEMARK_PATH=/home/dennu/Desktop/BRAKER/genemark-4.71 \
		--GUSHR_PATH=/home/dennu/Desktop/BRAKER/GUSHR \
		--PROTHINT_PATH=/home/dennu/Desktop/BRAKER/ProtHint/bin \
		--species braker_bathycoccus_prasinos \
		--softmasking
```

### SNAP training using maker

SNAP 2013-11-29 (https://github.com/KorfLab/SNAP) gene prediction model was trained through an first round of MAKER 2.31.9 (https://github.com/Yandell-Lab/maker) using assembled transcripts and plant protein database (plant_odb10) as evidences.

```bash
time mpiexec -n 12 maker -base maker_training_run {input.opts_ctl} {input.bopts_ctl} {input.exe_ctl}

maker2zff -x 0.25 -c 0.95 -d {input.maker_master_datastore}
fathom genome.ann genome.dna -gene-stats > gene-stats.log
fathom genome.ann genome.dna -validate > validate.log
fathom genome.ann genome.dna -categorize 1000 > categorize.log
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log

mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 
cd ..
hmm-assembler.pl genome params > snap_trained_modele.hmm
```

### Structural annotation

Merging of GENEMARK-ETP, AUGUSTUS and SNAP trained models was performed through a second round of MAKER 2.31.9 (https://github.com/Yandell-Lab/maker)

```bash
time mpiexec -n 12 maker -base maker_annotation_run {input.opts_ctl} {input.bopts_ctl} {input.exe_ctl}
```
