STRAIN_SUBSET = ["A8","C3"]

rule target_rule:
	input:
#		expand("/scratch/ldennu/merqury/{strain}_merqury/{strain}.qv", strain=STRAIN_SUBSET),
#		expand("/scratch/ldennu/merqury/{strain}_merqury/completeness.stats", strain=STRAIN_SUBSET),
		expand("/scratch/ldennu/busco/{strain}_busco/short_summary.specific.chlorophyta_odb10.{strain}_busco.txt", strain=STRAIN_SUBSET)

####################
rule merqury:
	input:
		merqury_singularity = "/data3/projects/bathycoccus/Tools_and_Scripts/Tools/Singularity.merqury.1.1.sif",
		illumina_reads1 = "/scratch/ldennu/illumina/{strain}_R1.fastq.gz",
		illumina_reads2 = "/scratch/ldennu/illumina/{strain}_R2.fastq.gz",
		assembly = "/scratch/ldennu/assembly/{strain}.fasta"		
	output:
		quality_stats = "/scratch/ldennu/merqury/{strain}_merqury/{strain}.qv",
		completeness_stats = "/scratch/ldennu/merqury/{strain}_merqury/completeness.stats"
	shell:
		"""
		module load system/singularity/3.6.0
		
		#mkdir /scratch/ldennu/merqury/{wildcards.strain}_merqury
		cd /scratch/ldennu/merqury/{wildcards.strain}_merqury

		singularity exec {input.merqury_singularity} meryl \
		count k=17 \
		{input.illumina_reads1} \
		output {wildcards.strain}_1.meryl

		singularity exec {input.merqury_singularity} meryl \
		count k=17 \
		{input.illumina_reads2} \
		output {wildcards.strain}_2.meryl

		singularity exec {input.merqury_singularity} meryl \
		union-sum {wildcards.strain}_1.meryl {wildcards.strain}_2.meryl \
		output {wildcards.strain}_allReads.meryl

		singularity exec {input.merqury_singularity} /opt/tools/merqury-1.1/merqury.sh \
		{wildcards.strain}_allReads.meryl/ \
		{input.assembly} \
		{wildcards.strain}
		"""
				
####################
rule busco:
	input:
		assembly = "/scratch/ldennu/assembly/{strain}.fasta",
		busco_singularity = "/data3/projects/bathycoccus/Tools_and_Scripts/Tools/busco_5.4.4.sif"
	output:
		busco_summary = "/scratch/ldennu/busco/{strain}_busco/short_summary.specific.chlorophyta_odb10.{strain}_busco.txt"
	shell:
		"""
		module load system/singularity/3.6.0

		cd /scratch/ldennu
		singularity exec {input.busco_singularity} busco \
		-i {input.assembly} \
		-l chlorophyta \
		-m genome \
		-c 12 \
		-o /busco/{wildcards.strain}_busco \
		-f
		"""

#singularity pull busco_5.4.4.sif docker://ezlabgva/busco:v5.4.4_cv1
####################
