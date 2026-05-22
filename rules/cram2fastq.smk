rule namesort_cram:
    input:
        cram = get_orig_cram_fh,
    output:
        cram = temp("data/intermediate_cram/{SAMPLE}.sorted.cram")
    params:
        tmpdir = "/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp/"
    resources:
        mem_mb = 48_000
    threads: 8
    shell:
        """
        module load samtools
        
        samtools sort -n \
                      -O CRAM \
                      -o {output.cram} \
                      -T {params.tmpdir} \
                      -@ {threads} \
                      -m 4G \
                      {input.cram}
        """


rule cram2fastq:
    input:
        cram = "data/intermediate_cram/{SAMPLE}.sorted.cram",
        ref = get_cram2fastq_ref,
    output:
        fq1 = temp("data/fastq/from_cram/{SAMPLE}.1.fastq.gz"),
        fq2 = temp("data/fastq/from_cram/{SAMPLE}.2.fastq.gz")
    threads: 8
    resources:
        mem_mb = 32_000
    shell:
        """
        module load samtools
        
        samtools fastq -@ {threads} \
                       --reference {input.ref} \
                       -1 {output.fq1} \
                       -2 {output.fq2} \
                       {input.cram}
        """