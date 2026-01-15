#!/bin/bash
set -e

fqs=(${snakemake_input[ora_list]})

for i in `seq 0 ${snakemake_params[n_input_files]}`; do
    echo $i
    fh="${fqs[$i]}"
    echo $fh
    ${snakemake_input[orad_binary]} --ora-reference /uufs/chpc.utah.edu/common/HIPAA/u1006375/src/orad.2.7.0.linux/oradata/ \
                                   -t ${snakemake[threads]} \
                                   -o data/fastq/${snakemake_wildcards[SAMPLE]}_${i}.fq.gz \
                                   --force \
                                   $fh
done
    

