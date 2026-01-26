#!/bin/bash
set -e

oras=(${snakemake_input[ora_list]})

for i in `seq 0 ${snakemake_params[n_input_files]}`; do
    fh="${oras[$i]}"
    ${snakemake_input[orad_binary]} --ora-reference /uufs/chpc.utah.edu/common/HIPAA/u1006375/src/orad.2.7.0.linux/oradata/ \
                                   -t ${snakemake[threads]} \
                                   -o ${snakemake_params[prefix]}data/fastq/${snakemake_wildcards[SAMPLE]}_${i}.fq.gz \
                                   --force \
                                   $fh
done
    

