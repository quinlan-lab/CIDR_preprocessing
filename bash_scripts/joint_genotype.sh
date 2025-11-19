#!/bin/bash
set -e

module load singularity

export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/deep_variant_tmp/

singularity exec --cleanenv -H $SINGULARITYENV_TMPDIR -B /usr/lib/locale/:/usr/lib/locale/ \
    ${snakemake_input[sif]} \
    /usr/local/bin/glnexus_cli \
    --dir ${snakemake_params[gl_nexus_prefix]} \
    --config DeepVariant_unfiltered \
    --mem-gbytes 64 \
    --threads ${snakemake[threads]} \
    ${snakemake_input[gvcfs]} > ${snakemake_output}