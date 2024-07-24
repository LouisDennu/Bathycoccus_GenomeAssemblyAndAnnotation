############################################

rule target_rule:
	input:
               	maker_gff = "/scratch/ldennu/maker_annotation_run.maker.output/maker_annotation_run.maker.gff",
                maker_transcripts_fasta = "/scratch/ldennu/maker_annotation_run.maker.output/maker_annotation_run.all.maker.transcripts.fasta",
                maker_protein_fasta = "/scratch/ldennu/maker_annotation_run.maker.output/maker_annotation_run.all.maker.proteins.fasta",
                maker_master_datastore = "/scratch/ldennu/maker_annotation_run.maker.output/maker_annotation_run_master_datastore_index.log"

rule maker_annotation_round:
	input:
		opts_ctl = "/scratch/ldennu/annotation_maker_opts.ctl",
		bopts_ctl = "/scratch/ldennu/maker_bopts.ctl",
		exe_ctl = "/scratch/ldennu/maker_exe.ctl",
		hmm_model = "/scratch/ldennu/snap/snap_trained_modele.hmm",
		GeneMark_ETP_model = "/scratch/ldennu/braker/GeneMark-ET/gmhmm.mod",
		Augustus_model = ""
	output:
		maker_gff = "/scratch/ldennu/maker_annotation_run.maker.output/maker_annotation_run.maker.gff",
		maker_transcripts_fasta = "/scratch/ldennu/maker_annotation_run.maker.output/maker_annotation_run.all.maker.transcripts.fasta",
		maker_protein_fasta = "/scratch/ldennu/maker_annotation_run.maker.output/maker_annotation_run.all.maker.proteins.fasta",
		maker_master_datastore = "/scratch/ldennu/maker_annotation_run.maker.output/maker_annotation_run_master_datastore_index.log"
	shell:
		"""
		module load system/python/3.8.12
		module load bioinfo/maker/2.31.9
		module load bioinfo/exonerate/2.4.0
		module load bioinfo/RepeatMasker/4.1.2
		module load bioinfo/blast/2.8.1+

		cd /scratch/ldennu
		time mpiexec -n 12 maker -base maker_annotation_run {input.opts_ctl} {input.bopts_ctl} {input.exe_ctl}
		cd /scratch/ldennu/maker_annotation_run.maker.output/
		gff3_merge -s -d ./maker_annotation_run_master_datastore_index.log > {output.maker_gff}
		fasta_merge -d ./maker_annotation_run_master_datastore_index.log
		"""

############################################
