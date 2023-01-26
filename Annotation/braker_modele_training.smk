rule target_rule:
	input:
		GeneMark_ET_model = "/home/dennu/Desktop/BRAKER/braker/GeneMark-ETP/gmhmm.mod",
		Augustus_model = "/home/dennu/Desktop/BRAKER/braker/species/braker_bathycoccus_prasinos"

rule braker_training:
	input:
		assembly = "/home/dennu/Desktop/BRAKER/4752.fasta.masked",
		bam = "/home/dennu/Desktop/BRAKER/4752_masked_RNAseq_HISAT2.sorted.bam",
		protein_db = "/home/dennu/Desktop/BRAKER/protein_db/plant_odb10.fasta"
	output:
		GeneMark_ET_model = "/home/dennu/Desktop/BRAKER/braker/GeneMark-ETP/gmhmm.mod",
		Augustus_model = directory("/home/dennu/Desktop/BRAKER/braker/species/braker_bathycoccus_prasinos")
	shell:
		"""

		export PATH=/home/dennu/Desktop/BRAKER/BRAKER/scripts/:$PATH
		export ETP=/home/dennu/Desktop/BRAKER/GeneMark-ETP/bin
		cd /home/dennu/Desktop/BRAKER	

		braker.pl \
		--genome {input.assembly} \
		--bam {input.bam} \
		--prot_seq {input.protein_db} \
		--etpmode \
		--cores 12 \
		--GENEMARK_PATH=/home/dennu/Desktop/BRAKER/genemark-4.71 \
		--GUSHR_PATH=/home/dennu/Desktop/BRAKER/GUSHR \
		--PROTHINT_PATH=/home/dennu/Desktop/BRAKER/ProtHint/bin \
		--species braker_bathycoccus_prasinos \
		--softmasking

		"""
