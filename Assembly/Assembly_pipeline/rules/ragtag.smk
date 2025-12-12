rule ragtag_scaffolding:
	input:
		pilon_fasta = lambda wc: f"{output_dir}/{wc.strain}/PILON/pilon.fasta",
		reference_genome = f"{reference_genome}"
	output:
		ragtag_folder = directory(f"{output_dir}/{{strain}}/RAGTAG"),
		ragtag_fasta = f"{output_dir}/{{strain}}/RAGTAG/ragtag.scaffold.fasta"
	log:
		out = f"{log_snakemake}/{{strain}}/out/ragtag.out",
		err = f"{log_snakemake}/{{strain}}/err/ragtag.err"
	threads: 1
	singularity: f"{singularity_dir}/ragtag.sif"
	conda: f"{singularity_dir}/env_yml/ragtag.yml"
	shadow: "shallow"
	shell:
		"""
		ragtag.py scaffold -Cr \
		-o {output.ragtag_folder} \
		{input.reference_genome} \
		{input.pilon_fasta} \
		1> {log.out} 2> {log.err}
		"""
