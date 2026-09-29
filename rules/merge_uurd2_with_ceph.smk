rule merge_fastq:
    input:
        old_fq = "data/fastq/from_cram/{SAMPLE}.{R}.fastq.gz",
        new_fq = "data/fastq/from_ora/{SAMPLE}.{R}.fastq.gz",
    output:
        fq = temp("data/fastq/merged/{SAMPLE}.{R}.fastq.gz"),
    shell:
        """
        cat {input.old_fq} {input.new_fq} > {output.fq}
        """