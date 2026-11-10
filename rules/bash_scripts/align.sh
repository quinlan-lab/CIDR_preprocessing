#!/bin/bash
set -euo pipefail

module load bwa/0.7.17 k8 gatk/4.6 samtools

SAMTOOLS_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp

T=${snakemake[threads]}
ALIGN_THREADS=$(( (T / 4) * 3 ))
SORT_THREADS=$(( T / 4 ))

bwa mem \
    -t $ALIGN_THREADS \
    -K ${snakemake_params[K]} \
    -R ${snakemake_params[read_group]} \
    -k 19 \
    ${snakemake_input[alignment_ref]} \
    ${snakemake_input[fq1]} \
    ${snakemake_input[fq2]} \
    | \
k8 ${snakemake_input[postalt]} \
    ${snakemake_input[alt_contigs]} | \
${snakemake_input[samblaster_binary]} \
    --addMateTags \
    --ignoreUnmated | \
samtools view -u - | \
samtools sort \
    -@ ${SORT_THREADS} \
    -o ${snakemake_output[cram]} \
    -m 4g \
    -O CRAM,version=3.0 \
    --output-fmt-option embed_ref=1 \
    --reference ${snakemake_input[cram_ref]} \
    --write-index \
    -T ${snakemake_params[tmp_prefix]} \
    -

