#!/bin/bash
set -e

module load singularity

export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/deep_variant_tmp/

singularity exec --cleanenv \
        -H $SINGULARITYENV_TMPDIR \
        -B /usr/lib/locale/:/usr/lib/locale/ \
            ${snakemake_input[sif]} \
            run_deepvariant \
                --model_type WGS \
                --num_shards ${snakemake[threads]} \
                --reads ${snakemake_input[cram]} \
                --output_vcf ${snakemake_output[vcf]} \
                --sample_name ${snakemake_wildcards[SAMPLE]} \
                --output_gvcf ${snakemake_output[gvcf]} \
                --ref ${snakemake_input[ref]} \
                --regions ${snakemake_wildcards[CHROM]} \
                ${snakemake_params[haploid_contigs_arg]}
                # --make_examples_extra_args "select_variant_types='snps',min_mapping_quality=1"