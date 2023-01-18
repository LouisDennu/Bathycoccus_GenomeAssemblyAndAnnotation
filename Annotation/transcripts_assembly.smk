STRAIN_SUBSET = ["4752"]

rule target_rule:
	input:
		expand ("/scratch/ldennu/RNAseq_snakemake/{strain}_masked_RNAseq_HISAT2_ST2_transcripts.gft", strain = STRAIN_SUBSET),
		expand ("/scratch/ldennu/RNAseq_snakemake/{strain}_masked_RNAseq_HISAT2_ST2_transcripts.fasta", strain = STRAIN_SUBSET),
		expand ("/scratch/ldennu/assembly/{strain}.fasta.masked", strain = STRAIN_SUBSET)

############################################
### MAKER Round 1

rule maker_round_1:
	input:
		opts.ctl = "/scratch/ldennu/training_maker_opts.ctl",
		bopts.ctl = "/scratch/ldennu/maker_bopts.ctl",
		exe.ctl = "/scratch/ldennu/maker_exe.ctl"
	output:
		maker_gff_round_1 = "/scratch/ldennu/
	shell:
		"""
		module load system/python/3.8.12
		module load bioinfo/maker/2.31.9
		module load bioinfo/exonerate/2.4.0
		module load bioinfo/RepeatMasker/4.1.2
		module load bioinfo/blast/2.8.1+

		cd /scratch/ldennu
		time mpiexec -n 12 maker -base Bprasinos_rnd1_011212022 round1_maker_opts.ctl maker_bopts.ctl maker_exe.ctl
		"""

############################################
### RNAseq to genome assembled transcripts

rule stringtie_transcripts_assembly:
	input:
		bam_file = "/scratch/ldennu/RNAseq_snakemake/{strain}_masked_RNAseq_HISAT2.sorted.bam",
		soft_masked_assembly = "/scratch/ldennu/assembly/{strain}.fasta.masked"
	output:
		transcripts_gft = "/scratch/ldennu/RNAseq_snakemake/{strain}_masked_RNAseq_HISAT2_ST2_transcripts.gft",
		transcripts_fasta = "/scratch/ldennu/RNAseq_snakemake/{strain}_masked_RNAseq_HISAT2_ST2_transcripts.fasta"
	shell:
		"""
		module load bioinfo/stringtie/1.3.4
		module load bioinfo/cufflinks/2.2.1-patched

		stringtie {input.bam_file} -p 12 -v -o {output.transcripts_gft}
		gffread {output.transcripts_gft} -g {input.soft_masked_assembly} -w {output.transcripts_fasta}
		"""

rule hisat2_mapping:
	input:
		index = "/scratch/ldennu/RNAseq_snakemake/{strain}_masked.index.8.ht2",
		RNAseq_illumina_1 = "/scratch/ldennu/RNAseq_snakemake/SRR1300453_1.fastq",
		RNAseq_illumina_2 = "/scratch/ldennu/RNAseq_snakemake/SRR1300453_2.fastq"
	output:
		sam_file = "/scratch/ldennu/RNAseq_snakemake/{strain}_masked_RNAseq_HISAT2.sam",
		bam_file = "/scratch/ldennu/RNAseq_snakemake/{strain}_masked_RNAseq_HISAT2.sorted.bam"
	shell:
		"""
		module load bioinfo/hisat2/2.1.0

		hisat2 -x /scratch/ldennu/RNAseq_snakemake/{wildcards.strain}_masked.index -1 {input.RNAseq_illumina_1} -2 {input.RNAseq_illumina_2} -S {output.sam_file} -p 12
		samtools view -@ 12 -Sb {output.sam_file} | samtools sort -o {output.bam_file}
		"""

rule hisat2_build:
	input:
		soft_masked_assembly = "/scratch/ldennu/assembly/{strain}.fasta.masked"
	output:
		index = "/scratch/ldennu/RNAseq_snakemake/{strain}_masked.index.8.ht2"
	shell:
		"""
		module load bioinfo/hisat2/2.1.0

		hisat2-build -p 12 {input.soft_masked_assembly} /scratch/ldennu/RNAseq_snakemake/{wildcards.strain}_masked.index
		"""

############################################
### Genome repeat masking

rule repeat_masker:
	input:
		assembly = "/scratch/ldennu/assembly/{strain}.fasta",
		repeat_db = "/scratch/ldennu/repeat_library/allNonRedundant9090BathyRepeatsFamilies.fasta"
	output:
		soft_masked_assembly = "/scratch/ldennu/assembly/{strain}.fasta.masked"
	shell:
		"""
		module load bioinfo/RepeatMasker/4.1.2

		cd /scratch/ldennu/assembly
		RepeatMasker -xsmall -lib {input.repeat_db} -s -pa 12 {input.assembly} --dir /scratch/ldennu/assembly
		"""

############################################
