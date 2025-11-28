assembly_pipeline_dir/
├── README.md
├── Snakefile
│
├── config/
│   ├── cluster.yaml
│   └── config.yaml
│
├── data/
│   ├── illumina/
│   │  ├── strainID_R1.fastq.gz
│   │  └── strainID_R2.fastq.gz
│   ├── ont/
│   │  └── strainID.fastq.gz
│   └── references/
│      └── reference.fasta/
│
├── logs/
│   └── snakemake/
│   └── cluster/
│
├── profiles/
│   └── slurm/
│      └── config.yaml
│
├── results/
│   └── StrainID/
│
├── rules/
│   ├── flye.smk
│   ├── medaka.smk
│   ├── pilon.smk
│   ├── ragtag.smk
│   └── tgs_gapcloser.smk
│
└── singularity/
    ├── env_yml/
    ├── flye.sif
    ├── medaka.sif
    ├── pilon.sif
    ├── ragtag.sif
    └── tgsgapcloser.sif
