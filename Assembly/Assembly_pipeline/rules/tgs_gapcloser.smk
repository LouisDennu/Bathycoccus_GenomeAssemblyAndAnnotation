rule tgs_gapcloser:
	input:
		ragtag_fasta = lambda wc: f"{output_dir}/{wc.strain}/RAGTAG/ragtag.scaffold.fasta",
		nanopore_reads = lambda wc: f"{ont_dir}/{wc.strain}.fastq.gz"
	output:
		nanopore_reads_fasta = f"{output_dir}/{{strain}}/TGS/{{strain}}_ONT.fasta",
		tgs_dir = directory(f"{output_dir}/{{strain}}/TGS"),
		tgs_gapcloser_assembly = f"{output_dir}/{{strain}}/TGS/{{strain}}.scaff_seqs"
	params:
		MODE = lambda wc: config["PILON_MODE"][wc.strain],
		r_round = config["r_round"],
		min_nread = config["min_nread"]
	log:
		out = f"{log_snakemake}/{{strain}}/out/tgs_gapcloser.out",
		err = f"{log_snakemake}/{{strain}}/err/tgs_gapcloser.err"
	threads: 4
	singularity: f"{singularity_dir}/tgsgapcloser.sif"
        conda: f"{singularity_dir}/env_yml/tgsgapcloser.yml"
	shadow: "full"
	shell:
		"""
		mkdir -p {output.tgs_dir}
		echo "[{wildcards.strain}] Converting Nanopore reads to FASTA..." 1> {log.out} 2> {log.err}
		seqtk seq -a {input.nanopore_reads} > {output.nanopore_reads_fasta} 2>> {log.err}
		
		echo "[{wildcards.strain}] Running TGS-GapCloser..." 1>> {log.out} 2>> {log.err}
		tgsgapcloser \
		--scaff {input.ragtag_fasta} \
		--reads {output.nanopore_reads_fasta} \
		--output {output.tgs_dir}/{wildcards.strain} \
		--samtools /opt/conda/envs/tgsgapcloser_env/bin/samtools \
		--racon /opt/conda/envs/tgsgapcloser_env/bin/racon \
		--r_round {params.r_round} \
		--p_round 0 \
		--min_nread {params.min_nread} \
		--thread {threads} \
		1>> {log.out} 2>> {log.err}
		"""
