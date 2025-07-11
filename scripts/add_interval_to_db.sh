#!/bin/bash
set -e

ls ${snakemake_params[db_dir]}

if [ -f "${snakemake_params[db_dir]}/callset.json" ]; then
  echo "GenomicsDB already created successfully"
else
  rmdir ${snakemake_params[db_dir]}
  echo "Removing the directory that Snakemake created automatically in order to initialize a new GenomicsDB"
fi

module load gatk/4.6
        
gatk --java-options '-Xmx32g -Xms32g' \
    GenomicsDBImport \
    ${snakemake_params[workspace_arg]} ${snakemake_params[db_dir]} \
    --batch-size 24 \
    --sample-name-map ${snakemake_input[sample_map]} \
    --tmp-dir /scratch/ucgd/lustre-labs/quinlan/u1006375/gatk_tmp \
    --reader-threads ${snakemake[threads]} \
    ${snakemake_params[interval_arg]}

echo "done" >> ${snakemake_output[completed]}
