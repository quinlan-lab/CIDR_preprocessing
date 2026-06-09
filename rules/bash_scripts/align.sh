#!/bin/bash
set -euo pipefail

module load gatk/4.6

export SAMTOOLS_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp

${snakemake_input[bwa_binary]} mem \
    -t ${snakemake[threads]} \
    -K ${snakemake_params[K]} \
    -R ${snakemake_params[read_group]} \
    ${snakemake_input[ref]} \
    ${snakemake_input[fq1]} \
    ${snakemake_input[fq2]} \
    | \
  gatk SortSam \
    --java-options -Xmx32g \
    --MAX_RECORDS_IN_RAM ${snakemake_params[max_records_in_ram]} \
    -I /dev/stdin \
    -O ${snakemake_output[bam]} \
    --SORT_ORDER coordinate \
    --TMP_DIR $SAMTOOLS_TMPDIR

