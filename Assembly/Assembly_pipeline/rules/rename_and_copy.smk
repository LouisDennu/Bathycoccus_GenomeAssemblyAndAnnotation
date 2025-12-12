rule rename_and_copy:
	input:
		TGS_assembly_file = f"{output_dir}/{{strain}}/TGS/{{strain}}.scaff_seqs"
	output:
		final_assembly = f"{output_dir}/{{strain}}/final/{{strain}}_assembly.fasta"
	threads: 1
	shell:
		"""
		mkdir -p $(dirname {output.final_assembly})
		rsync -vaurL {input.TGS_assembly_file} {output.final_assembly}
		"""
