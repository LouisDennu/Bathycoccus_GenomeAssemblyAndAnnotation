rule flye_draft_assembly:
	input:
		nanopore_reads = lambda wc: f"{ont_dir}/{wc.strain}.fastq.gz"
	output:
		flye_folder = directory(f"{output_dir}/{{strain}}/FLYE"),
		flye_assembly = f"{output_dir}/{{strain}}/FLYE/assembly.fasta"
	params:
		genome_size = config["genome_size"]
	log:
		out = f"{log_snakemake}/{{strain}}/out/flye.out",
		err = f"{log_snakemake}/{{strain}}/err/flye.err"
	threads: 8
	singularity: f"{singularity_dir}/flye.sif"
	conda: f"{singularity_dir}/env_yml/flye.yml"
	shadow: "shallow"
	shell:
		"""
		flye --nano-hq {input.nanopore_reads} \
		--genome-size {params.genome_size} \
		--scaffold \
		--threads {threads} \
		--out-dir {output.flye_folder} \
		1> {log.out} 2> {log.err}
		"""
