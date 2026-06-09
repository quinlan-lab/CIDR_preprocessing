# rule namesort_cram:
#     input:
#         cram = get_orig_cram_fh,
#     output:
#         cram = temp("data/intermediate_cram/{SAMPLE}.sorted.cram")
#     params:
#         tmpdir = "/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp/",
#         mem = "4G"
#     resources:
#         mem_mb = 72_000,
#         runtime = 1440,
#     threads: 12
#     shell:
#         """
#         module load samtools
        
#         samtools sort -n \
#                       -O CRAM \
#                       -o {output.cram} \
#                       -T {params.tmpdir} \
#                       -@ {threads} \
#                       -m {params.mem} \
#                       {input.cram}
#         """


# rule cram2fastq:
#     input:
#         cram = "data/intermediate_cram/{SAMPLE}.sorted.cram",
#         ref = get_cram2fastq_ref,
#     output:
#         fq1 = temp("data/fastq/from_cram/{SAMPLE}.1.fastq.gz"),
#         fq2 = temp("data/fastq/from_cram/{SAMPLE}.2.fastq.gz")
#     threads: 16
#     resources:
#         mem_mb = 32_000
#     shell:
#         """
#         module load samtools
        
#         samtools fastq -@ {threads} \
#                        --reference {input.ref} \
#                        -1 {output.fq1} \
#                        -2 {output.fq2} \
#                        {input.cram}
#         """

# https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
rule collate_fastq:
    input:
        cram = get_orig_cram_fh,
        ref = get_cram2fastq_ref,
        samtools_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/samtools"
    output:
        fq1 = temp("data/fastq/from_cram/{SAMPLE}.1.fastq.gz"),
        fq2 = temp("data/fastq/from_cram/{SAMPLE}.2.fastq.gz"),
        fq0 = temp("data/fastq/from_cram/{SAMPLE}.0.fastq.gz")
    threads: 16
    resources:
        mem_mb = 32_000,
        runtime = 1440
    params:
        tmpdir = "/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp/",

    shell:
        """        
        export TMPDIR={params.tmpdir}
        echo $TMPDIR
        {input.samtools_binary} collate -@ {threads} \
                         --reference {input.ref} \
                         -Oun64 {input.cram} \
            | samtools fastq -@ {threads} -0 {output.fq0} -1 {output.fq1} -2 {output.fq2} -
        """
