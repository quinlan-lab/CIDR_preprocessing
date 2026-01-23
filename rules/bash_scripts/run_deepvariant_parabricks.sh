#!/bin/bash
set -e

pwd
nvidia-smi

module load singularity

export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/parabricks_tmp/

singularity exec --cleanenv \
        --nv \
        -H $SINGULARITYENV_TMPDIR \
        --bind $(pwd):/workdir \
        --bind $(pwd):/outputdir \
        --pwd /workdir \
        ${snakemake_input[sif]} \
        pbrun deepvariant \
            --gvcf \
            --in-bam ${snakemake_input[cram]} \
            --out-variants ${snakemake_output[gvcf]} \
            --ref ${snakemake_input[ref]} \
            --tmp-dir $SINGULARITYENV_TMPDIR \
            ${snakemake_params[haploid_contigs_arg]}
