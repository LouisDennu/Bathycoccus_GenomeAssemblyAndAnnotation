# ONT reads basecalling

## Basecalling

Basecalling of ONT data was performed with Guppy 6.1.2 (https://nanoporetech.com) using model dna_r9.4.1_450bps_sup.

```bash
guppy_basecaller -c dna_r9.4.1_450bps_sup.cfg -i ${INPUT} -r -s ${OUTPUT} --barcode_kits "SQK-RBK004" --compress_fastq --num_callers 8 --gpu_runners_per_device 8 --min_qscore 7 -x cuda:${CUDA}
```

## Quality control

Quality control was performed using NanoPlot 1.19.0 (https://github.com/wdecoster/NanoPlot) to obtain statistics with graphical representation for basecalled data.

```bash
NanoPlot -t 4 -o ./nanoplot --summary ./sequencing_summary.txt
```
