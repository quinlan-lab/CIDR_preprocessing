#!/bin/bash
set -euo pipefail

module load bwa/0.7.19 gatk/4.6

bwa mem \
    -t ${snakemake[threads]} \
    -K 96000000 \
    -R '@RG\tID:{wildcards.SAMPLE}\tLB:lib1\tPL:ILLUMINA\tSM:{wildcards.SAMPLE}\tPU:{wildcards.SAMPLE}' \
    ${snakemake_input[ref]} \
    ${snakemake_input[fq1]} \
    ${snakemake_input[fq2]} \
    | \
  gatk SortSam \
    --java-options -Xmx32g \
    --MAX_RECORDS_IN_RAM 5000000 \
    -I /dev/stdin \
    -O ${snakemake_output[bam]} \
    --SORT_ORDER coordinate

