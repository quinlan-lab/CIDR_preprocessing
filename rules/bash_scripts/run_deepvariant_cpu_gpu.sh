#!/bin/bash
set -e


gpu=${snakemake_params[gpu]}
echo $gpu
module load singularity

if [ "$gpu" -eq "1" ]; then

    export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/parabricks_tmp/

    singularity exec --cleanenv \
            --nv \
            -H $SINGULARITYENV_TMPDIR \
            --bind $(pwd):/workdir \
            --bind $(pwd):/outputdir \
            --pwd /workdir \
            ${snakemake_params[sif]} \
            pbrun deepvariant \
                --gvcf \
                --in-bam ${snakemake_input[cram]} \
                --out-variants ${snakemake_output[gvcf]} \
                --ref ${snakemake_input[ref]} \
                --tmp-dir $SINGULARITYENV_TMPDIR \
                -L ${snakemake_wildcards[CHROM]} \
                --num-cpu-threads-per-stream 2 \
                ${snakemake_params[haploid_contigs_arg]}
else

    export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/deep_variant_tmp_2/

    singularity exec --cleanenv \
            -H $SINGULARITYENV_TMPDIR \
            -B /usr/lib/locale/:/usr/lib/locale/ \
                ${snakemake_params[sif]} \
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
fi