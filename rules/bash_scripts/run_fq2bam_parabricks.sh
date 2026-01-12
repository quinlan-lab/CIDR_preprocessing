#!/bin/bash
set -e

pwd
nvidia-smi

module load singularity

export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/parabricks_tmp/
export SAMTOOLS_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp

singularity exec --cleanenv \
        --nv \
        -H $SINGULARITYENV_TMPDIR \
        --bind $(pwd):/workdir \
        --bind $(pwd):/outputdir \
        --pwd /workdir \
        ${snakemake_input[sif]} \
        pbrun fq2bam \
            --ref ${snakemake_input[ref]} \
            --in-fq ${snakemake_input[fq1]} ${snakemake_input[fq2]} \
            --out-bam ${snakemake_output[bam]} \
            --read-group-sm ${snakemake_wildcards[SAMPLE]} \
            --read-group-lb lib1 \
            --read-group-pl ILLUMINA \
            --read-group-id-prefix unit1 \
            --memory-limit 48 \
            --tmp-dir $SAMTOOLS_TMPDIR \
            --bwa-options="-K 10000000"