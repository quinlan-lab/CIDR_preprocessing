#!/bin/bash
set -e

n=${snakemake_params[n_input_files]}
echo ${snakemake_params[n_input_files]}
echo $n
echo ${snakemake_params[fastq1_list]}
if [ "$n" -eq "1" ]; then
    echo "yes"
    mv ${snakemake_params[fastq1_list]:0} ${snakemake_output[fq1]}
    mv ${snakemake_params[fastq2_list]:0} ${snakemake_output[fq2]}
else
    cat ${snakemake_params[fastq1_list]} > ${snakemake_output[fq1]}
    cat ${snakemake_params[fastq2_list]} > ${snakemake_output[fq2]}

    rm ${snakemake_params[fastq1_list]}
    rm ${snakemake_params[fastq2_list]}
fi
