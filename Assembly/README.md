# Bathycoccus prasinos Pangenome
## Genome de novo assembly and polishing

![alt text](https://github.com/LouisDennu/BathycoccusPangenome/blob/main/Assembly/Bathycoccus_pangenome_pipelines-Genome_assembly_pipeline.png)

---

### Genome assembly

Genome assembly was produced with ONT reads by FLYE 2.9.

https://github.com/fenderglass/Flye

```shell
python flye --nano-hq ${FASTQ_ONT} --genome-size 15m --scaffold --threads 6 --out-dir ${OUTDIR}
```

### Assembly polishing

Assembly polishing using nanopore reads was done using MEDAKA 1.5.

https://github.com/nanoporetech/medaka

```shell
medaka_consensus -i ${FASTQ_ONT} -d ${ASSEMBLY} -o ${OUTDIR} -m r941_min_high_g360 -t 2
```

### Assembly correction

Assembly correction using illumina reads was done using BWA 0.7.17 for read mapping, Samtools 1.9 and PILON 1.24.

https://github.com/lh3/bwa
https://github.com/samtools/samtools
https://github.com/broadinstitute/pilon

```shell
bwa index ${ASSEMBLY}
bwa mem ${ASSEMBLY} ${FASTQ_ILLUMINA_1} ${FASTQ_ILLUMINA2} | samtools view -b | samtools sort -o ${ILLUMINA.sorted.bam}
samtools index ${ILLUMINA.sorted.bam}
java -Xmx16G -jar pilon-1.24.jar --genome ${ASSEMBLY}	--frags	${ILLUMINA.sorted.bam} --outdir ${OUTDIR}
```

### Scaffolding and direction of contigs on the reference genome

Scaffolding and direction of contigs on the reference genome (Strain RCC1105) was done using RagTag 2.1.0.
Unmapped contigs were grouped in Chr_0, considered as contaminations and discarded.

https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_002220235.1/
https://github.com/malonge/RagTag

```shell
ragtag.py scaffold -Cr ${REFERENCE_GENOME} ${ASSEMBLY}
```

### Gap filling of scaffolded assemblt

Gap filling was done using TGS-GapCloser 1.0.1 with Samtools 1.9 and PILON 1.24. Seqtk 1.3-r106 is used to convert FASTQ reads to FASTA for TGS-GapCloser input.
ONT reads are corrected by Illumina reads using PILON and mapped on the assembly, gap crossing reads are used to fill the gap.
One round of Pilon correction is used to stay consistent with global assembly correction.

https://github.com/lh3/seqtk
https://github.com/samtools/samtools
https://github.com/broadinstitute/pilon
https://github.com/BGI-Qingdao/TGS-GapCloser

```shell
seqtk seq -a ${FASTQ_ONT} > ${FASTA_ONT}
TGS-GapCloser.sh --scaff ${ASSEMBLY} --reads ${FASTA_ONT} --output ${OUTDIR} --pilon pilon-1.24.jar --java /usr/bin/java --ngs ${FASTQ_ILLUMINA_1} ${FASTQ_ILLUMINA2} --samtools /usr/local/samtools-1.9/bin/samtools --p_round 1 --r_round 0 --min_nread 3 --pilon_mem 16G
```

---

## Quality control of assembly

### Genome statistics

Genome statistics were produced using Assembly-stats 1.0.1.
https://github.com/sanger-pathogens/assembly-stats

```shell
assembly-stats -t ${ASSEMBLY} > ${STATS_ASSEMBLY}
```

### Genome completion and quality assesment

Genome completion assesment through search of conserved core genes was done using BUSCO 5.4.4 singularity image.

https://busco.ezlab.org/

Genome completion and quality assesment through Illumina reads mapping was done using MERQURY 1.1 singularity image.

https://github.com/marbl/merqury

```shell
singularity exec ${BUSCO_SINGULARITY} busco -i ${ASSEMBLY} -l chlorophyta -m genome -c 12 -o ${OUTDIR} -f
```

```shell
singularity exec ${MERQURY_SINGULARITY} meryl \
		count k=17 \
		${FASTQ_ILLUMINA_1} \
		output ${ILLUMINA_1.MERYL}

singularity exec ${MERQURY_SINGULARITY} meryl \
		count k=17 \
		${FASTQ_ILLUMINA_2} \
		output ${ILLUMINA_2.MERYL}

singularity exec ${MERQURY_SINGULARITY} meryl \
		union-sum ${ILLUMINA_1.MERYL} ${ILLUMINA_2.MERYL} \
		output ${ILLUMINA_ALLREADS.MERYL}

singularity exec ${MERQURY_SINGULARITY} /opt/tools/merqury-1.1/merqury.sh \
		${ILLUMINA_ALLREADS.MERYL} \
		${ASSEMBLY} \
		${STRAIN_NAME}
```

---
