rule target_rule:
	input:
		GeneMark_ET_model = "/scratch/ldennu/braker/GeneMark-ET/gmhmm.mod"

rule braker_training:
	input:
		braker_singularity = "/scratch/ldennu/braker_singularity/braker3.sif",
		assembly = "/scratch/ldennu/assembly/4752.fasta.masked",
		bam = "/scratch/ldennu/RNAseq_snakemake/4752_masked_RNAseq_HISAT2.sorted.bam",
	output:
		GeneMark_ET_model = "/scratch/ldennu/braker/GeneMark-ET/gmhmm.mod"
	shell:
		"""
		module load system/singularity/3.6.0

		cd /scratch/ldennu/	

		singularity exec {input.braker_singularity} braker.pl \
		--genome {input.assembly} \
		--bam {input.bam} \
		--threads 12 \
		--GENEMARK_PATH=/home/dennu/tools/genemark-4.71 \
 		--softmasking \
 		--UTR on
		"""
