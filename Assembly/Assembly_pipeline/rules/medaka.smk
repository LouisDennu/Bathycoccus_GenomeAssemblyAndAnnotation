rule medaka_correction:
	input:
		nanopore_reads = lambda wc: f"{ont_dir}/{wc.strain}.fastq.gz",
		flye_assembly = lambda wc: f"{output_dir}/{wc.strain}/FLYE/assembly.fasta"
	output:
		medaka_folder = directory(f"{output_dir}/{{strain}}/MEDAKA"),
		medaka_consensus = f"{output_dir}/{{strain}}/MEDAKA/consensus.fasta"
	params:
		medaka_model = config["medaka_model"]
	log:
		out = f"{log_snakemake}/{{strain}}/out/medaka.out",
		err = f"{log_snakemake}/{{strain}}/err/medaka.err"
	threads: 4
	singularity: f"{singularity_dir}/medaka.sif"
        conda: f"{singularity_dir}/env_yml/medaka.yml"
	shadow: "full"
	shell:
		"""
		medaka_consensus \
		-i {input.nanopore_reads} \
		-d {input.flye_assembly} \
		-o {output.medaka_folder} \
		-m {params.medaka_model} \
		-t {threads} \
		1> {log.out} 2> {log.err}
		"""
