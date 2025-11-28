How to setup:

 - Modify base_dir and log_dir variable to the actual PATH of the pipeline directory in ./config/config.yaml
 - Build singularity image via ./build_singularity.sh if able
 - Put ONT read files in ./data/ont
 - Put optional illumina read files in ./data/illumina
 - Modify reference genome for scaffolding in ./data/references if necessary

How to run:
 - For local use, run the pipeline with command "snakemake -j X -p --use-singularity" or "snakemake -j X -p --use-conda" with X the number of Job
 - For HPC use with SLURM scheduler, run submit_slurm.sh (modify partition names, module names and other parameter accordingly)

Assembly_pipeline/
- logs
- config
  - config.yaml
  - cluster.yaml
- README.md
- Snakefile
- build_singularity.sh
- submit_slurm.sh
- profiles
  - slurm
    - config.yaml
- rules
  - ragtag.smk
  - flye.smk
  - pilon.smk
  - medaka.smk
  - tgs_gapcloser.smk
- results
- singularity
  - env_yml
    - medaka.yml
    - ragtag.yml
    - tgsgapcloser.yml
    - flye.yml
    - pilon.yml
  - singu.def
- data
  - ont
  - references
    - RCC4222_assembly.fasta
  - illumina
