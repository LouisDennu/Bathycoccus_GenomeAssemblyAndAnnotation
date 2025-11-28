rule pilon_polishing:
	input:
                illumina_r1 = lambda wc: config["ILLUMINA_READS"][wc.strain][0] if config["PILON_MODE"][wc.strain] == "polish" else [],
                illumina_r2 = lambda wc: config["ILLUMINA_READS"][wc.strain][1] if config["PILON_MODE"][wc.strain] == "polish" else [],
		medaka_consensus = lambda wc: f"{output_dir}/{wc.strain}/MEDAKA/consensus.fasta"
	output:
		bam_file = f"{output_dir}/{{strain}}/{{strain}}_FLYE_MEDAKA.sorted.bam",
		pilon_fasta = f"{output_dir}/{{strain}}/PILON/pilon.fasta",
		pilon_folder = directory(f"{output_dir}/{{strain}}/PILON")
	params:
		MODE = lambda wc: config["PILON_MODE"][wc.strain],
		pilon_memory = config["pilon_memory"]
	log:
		out = f"{log_snakemake}/{{strain}}/out/pilon.out",
		err = f"{log_snakemake}/{{strain}}/err/pilon.err"
	threads: 4
	singularity: f"{singularity_dir}/pilon.sif"
        conda: f"{singularity_dir}/env_yml/pilon.yml"
	shadow: "full"
	shell:
		"""
		MODE="{params.MODE}"

		mkdir -p {output.pilon_folder}

		if [ "$MODE" = "skip" ]; then
			echo "[PILON] Skipping polishing, copying medaka consensus" 1> {log.out} 2> {log.err}
			cp {input.medaka_consensus} {output.pilon_fasta}
			touch {output.bam_file}
		else

			echo "[PILON] {wildcards.strain}: POLISHING with Illumina" 1> {log.out} 2> {log.err}

			bwa index {input.medaka_consensus} 1>> {log.out} 2>> {log.err}

			(
			bwa mem -t {threads} {input.medaka_consensus} {input.illumina_r1} {input.illumina_r2} |
			samtools view -b |
			samtools sort -o {output.bam_file}
			) 1>> {log.out} 2>> {log.err}

			samtools index {output.bam_file} 1>> {log.out} 2>> {log.err}

			pilon \
				--genome {input.medaka_consensus} \
				--frags {output.bam_file} \
				--outdir {output.pilon_folder} \
				-Xmx{params.pilon_memory} \
				1>> {log.out} 2>> {log.err}
		fi
		"""


