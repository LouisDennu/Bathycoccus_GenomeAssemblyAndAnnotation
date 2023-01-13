STRAIN_SUBSET = ["A1", "A7", "C218", "C3", "G8"]

rule target_rule:
	input:
		expand("/scratch/ldennu/Bathycoccus_assembly/{strain}/TGS/GC.scaff_seqs", strain = STRAIN_SUBSET)

rule tgs_gapcloser:
        input:
                tgs_gapcloser = "/home/dennu/tools/TGS-GapCloser/TGS-GapCloser.sh",
                ragtag_folder = "/scratch/ldennu/Bathycoccus_assembly/{strain}/RAGTAG",
                nanopore_reads = "/scratch/ldennu/ONT/{strain}.fastq.gz",
                illumina_reads1 = "/scratch/ldennu/illumina/{strain}_R1.fastq.gz",
                illumina_reads2 = "/scratch/ldennu/illumina/{strain}_R2.fastq.gz",
                pilon = "/home/dennu/tools/pilon-1.24.jar"
        output:
                nanopore_reads_fasta = "/scratch/ldennu/Bathycoccus_assembly/{strain}/TGS/{strain}_ONT.fasta",
                tgs_gapcloser_assembly = "/scratch/ldennu/Bathycoccus_assembly/{strain}/TGS/GC.scaff_seqs"
        shell:
                """
                module load bioinfo/seqtk/1.3-r106

                mkdir -p /scratch/ldennu/Bathycoccus_assembly/{wildcards.strain}/TGS
                cd /scratch/ldennu/Bathycoccus_assembly/{wildcards.strain}/TGS
                seqtk seq -a {input.nanopore_reads} > {output.nanopore_reads_fasta}
                
                {input.tgs_gapcloser} --scaff {input.ragtag_folder}/ragtag.scaffold.fasta \
                --reads {output.nanopore_reads_fasta} \
                --output GC \
                --pilon {input.pilon} \
                --java /usr/java/default/bin/java \
                --ngs {input.illumina_reads1} {input.illumina_reads2} \
                --samtools /usr/local/samtools-1.9/bin/samtools \
                --p_round 1 \
                --r_round 0 \
                --min_nread 3 \
                --pilon_mem 16G
                """

rule ragtag_scaffolding:
        input:
                pilon_folder = "/scratch/ldennu/Bathycoccus_assembly/{strain}/PILON",
                reference_genome = "/home/dennu/data/ASSEMBLY/4752_Assembly_files/4752_RB8_FLYE_MEDAKA_PILON_RAGTAG_GC.fasta"
        output:
                ragtag_folder = directory("/scratch/ldennu/Bathycoccus_assembly/{strain}/RAGTAG")
        shell:
                """
                module load bioinfo/ragtag/2.1.0
                ragtag.py scaffold -Cr -o {output.ragtag_folder} {input.reference_genome} {input.pilon_folder}/pilon.fasta
                """

rule pilon_polishing:
        input:
                pilon = "/home/dennu/tools/pilon-1.24.jar",
                illumina_reads1 = "/scratch/ldennu/illumina/{strain}_R1.fastq.gz",
                illumina_reads2 = "/scratch/ldennu/illumina/{strain}_R2.fastq.gz",
                medaka_folder = "/scratch/ldennu/Bathycoccus_assembly/{strain}/MEDAKA"
        output:
                bam_file = "/scratch/ldennu/Bathycoccus_assembly/{strain}/{strain}_FLYE_MEDAKA.sorted.bam",
                pilon_folder = directory("/scratch/ldennu/Bathycoccus_assembly/{strain}/PILON")
        shell:
                """
                module load bioinfo/bwa/0.7.17
                module load bioinfo/samtools/1.9
                
                bwa index {input.medaka_folder}/consensus.fasta
                bwa mem {input.medaka_folder}/consensus.fasta {input.illumina_reads1} {input.illumina_reads2} | samtools view -b | samtools sort -o {output.bam_file}
                samtools index {output.bam_file}
                java -Xmx16G -jar {input.pilon} --genome {input.medaka_folder}/consensus.fasta --frags {output.bam_file} --outdir {output.pilon_folder}
                """

rule medaka_correction:
        input:
                nanopore_reads = "/scratch/ldennu/ONT/{strain}.fastq.gz",
                flye_folder = "/scratch/ldennu/Bathycoccus_assembly/{strain}/FLYE"
        output:
                medaka_folder = directory("/scratch/ldennu/Bathycoccus_assembly/{strain}/MEDAKA")
        shell:
                """
                module load bioinfo/medaka/1.5
                medaka_consensus -i {input.nanopore_reads} -d {input.flye_folder}/assembly.fasta -o {output.medaka_folder} -m r941_min_high_g360 -t 2
                """

rule flye_draft_assembly:
        input:
              	flye = "/home/dennu/tools/Flye/bin/flye",
                nanopore_reads = "/scratch/ldennu/ONT/{strain}.fastq.gz"
        output:
               	flye_folder = directory("/scratch/ldennu/Bathycoccus_assembly/{strain}/FLYE")
        shell:
              	"""
                python {input.flye} --nano-hq {input.nanopore_reads} --genome-size 15m --scaffold --threads 6 --out-dir {output.flye_folder}
                """

