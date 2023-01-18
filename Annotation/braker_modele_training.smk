STRAIN_SUBSET = [""]

rule target_rule:
	input:
		

rule braker_training:
	input:
		assembly = ""
		bam = ""
	output:
	
	shell:
		"""
		module load system/Miniconda3/1.0
		module load system/perl/5.24.0
		module unload system/perl/5.24.0 
		module load bioinfo/BRAKER/2.1.5

		conda activate genemark-ET_env

		perl /home/dennu/tools/BRAKER-2.1.6/scripts/braker.pl --genome {input.assembly} --bam {input.bam} --cores 12 --GENEMARK_PATH /home/dennu/tools/genemark --AUGUSTUS_CONFIG_PATH /home/dennu/tools/augustus-3.3.3/config
		"""
