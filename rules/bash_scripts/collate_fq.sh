#!/bin/bash
set -euo pipefail

module load samtools

samtools collate \
         -@ ${snakemake[threads]} \
         --reference ${snakemake_input[ref]} \
         -O \
         -u \
         -n32 \
         ${snakemake_input[cram]} \
         ${snakemake_params[tmp_prefix]} | \
samtools fastq -@ ${snakemake[threads]} \
               -0 ${snakemake_output[fq0]} \
               -1 ${snakemake_output[fq1]} \
               -2 ${snakemake_output[fq2]} -