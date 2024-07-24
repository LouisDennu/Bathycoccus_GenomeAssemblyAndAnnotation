############################################
### MAKER Round 1

rule target_rule:
	input:
               	maker_training_gff = "/scratch/ldennu/maker_training_run.maker.output/maker_training_run.maker.gff",
                maker_training_transcripts_fasta = "/scratch/ldennu/maker_training_run.maker.output/maker_training_run.all.maker.transcripts.fasta",
                maker_training_protein_fasta = "/scratch/ldennu/maker_training_run.maker.output/maker_training_run.all.maker.proteins.fasta",
		hmm = "/scratch/ldennu/snap/snap_trained_modele.hmm"

rule SNAP_training:
	input:
		maker_master_datastore = "/scratch/ldennu/maker_training_run.maker.output/maker_training_run_master_datastore_index.log"
	output:
		hmm = "/scratch/ldennu/snap/snap_trained_modele.hmm"
	shell:
		"""
		module load bioinfo/snap/2013-11-29
		
		cd /scratch/ldennu/snap
		maker2zff -x 0.25 -c 0.95 -d {input.maker_master_datastore}
		fathom genome.ann genome.dna -gene-stats > gene-stats.log
		fathom genome.ann genome.dna -validate > validate.log
		fathom genome.ann genome.dna -categorize 1000 > categorize.log
		fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log

		mkdir params
		cd params
		forge ../export.ann ../export.dna > ../forge.log 
		cd ..
		hmm-assembler.pl genome params > snap_trained_modele.hmm
		"""

rule maker_training_round:
	input:
		opts_ctl = "/scratch/ldennu/training_maker_opts.ctl",
		bopts_ctl = "/scratch/ldennu/maker_bopts.ctl",
		exe_ctl = "/scratch/ldennu/maker_exe.ctl"
	output:
		maker_training_gff = "/scratch/ldennu/maker_training_run.maker.output/maker_training_run.maker.gff",
		maker_training_transcripts_fasta = "/scratch/ldennu/maker_training_run.maker.output/maker_training_run.all.maker.transcripts.fasta",
		maker_training_protein_fasta = "/scratch/ldennu/maker_training_run.maker.output/maker_training_run.all.maker.proteins.fasta",
		maker_master_datastore = "/scratch/ldennu/maker_training_run.maker.output/maker_training_run_master_datastore_index.log"
	shell:
		"""
		module load system/python/3.8.12
		module load bioinfo/maker/2.31.9
		module load bioinfo/exonerate/2.4.0
		module load bioinfo/RepeatMasker/4.1.2
		module load bioinfo/blast/2.8.1+

		cd /scratch/ldennu
		time mpiexec -n 12 maker -base maker_training_run {input.opts_ctl} {input.bopts_ctl} {input.exe_ctl}
		cd /scratch/ldennu/maker_training_run.maker.output/
		gff3_merge -s -d ./maker_training_run_master_datastore_index.log > {output.maker_training_gff}
		fasta_merge -d ./maker_training_run_master_datastore_index.log
		"""

############################################
