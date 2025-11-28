rule tgs_gapcloser:
	input:
		ragtag_fasta = lambda wc: f"{output_dir}/{wc.strain}/RAGTAG/ragtag.scaffold.fasta",
		nanopore_reads = lambda wc: f"{ont_dir}/{wc.strain}.fastq.gz",
		illumina_r1 = lambda wc: config["ILLUMINA_READS"][wc.strain][0] if config["PILON_MODE"][wc.strain] == "polish" else [],
		illumina_r2 = lambda wc: config["ILLUMINA_READS"][wc.strain][1] if config["PILON_MODE"][wc.strain] == "polish" else []
	output:
		nanopore_reads_fasta = f"{output_dir}/{{strain}}/TGS/{{strain}}_ONT.fasta",
		tgs_dir = directory(f"{output_dir}/{{strain}}/TGS"),
		tgs_gapcloser_assembly = f"{output_dir}/{{strain}}/TGS/{{strain}}.scaff_seqs"
	params:
		MODE = lambda wc: config["PILON_MODE"][wc.strain],
		pilon_memory = config["pilon_memory"],
		r_round = config["r_round"],
		p_round = config["p_round"],
		min_nread = config["min_nread"]
	log:
#		out = f"{log_snakemake}/{{strain}}/out/tgs_gapcloser.out",
		err = f"{log_snakemake}/{{strain}}/err/tgs_gapcloser.err"
	threads: 4
	singularity: f"{singularity_dir}/tgsgapcloser.sif"
	shadow: "shallow"
	shell:
		"""
		MODE="{params.MODE}"
		mkdir -p {output.tgs_dir}
		seqtk seq -a {input.nanopore_reads} > {output.nanopore_reads_fasta}

	        if [ "$MODE" = "skip" ]; then

			echo "[TGS] {wildcards.strain}: SKIP mode → running without Illumina polishing"
			
			tgsgapcloser \
			--scaff {input.ragtag_fasta} \
			--reads {output.nanopore_reads_fasta} \
			--output {output.tgs_dir}/{wildcards.strain} \
			--samtools $(which samtools) \
			--racon $(which racon) \
			--r_round {params.r_round} \
			--p_round 0 \
			--min_nread {params.min_nread} \
			--thread {threads}
		else

			echo "[TGS] {wildcards.strain}: POLISH mode → running with Illumina polishing"
			PILON_JAR=$(find "$(dirname "$(dirname "$(which pilon)")")" -type f -name "pilon*.jar" | head -n 1)

			tgsgapcloser \
			--scaff {input.ragtag_fasta} \
			--reads {output.nanopore_reads_fasta} \
			--output {output.tgs_dir}/{wildcards.strain} \
			--pilon $PILON_JAR\
			--java $(which java) \
			--ngs {input.illumina_r1} {input.illumina_r2} \
			--samtools $(which samtools) \
			--p_round {params.p_round} \
			--r_round 0 \
			--min_nread {params.min_nread} \
			--pilon_mem 32G \
			--thread {threads}
		fi
		"""
