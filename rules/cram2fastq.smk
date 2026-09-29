# https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
rule collate_fastq:
    input:
        cram = get_orig_cram_fh,
        ref = get_cram2fastq_ref
    output:
        fq1 = temp("data/fastq/from_cram/{SAMPLE}.1.fastq.gz"),
        fq2 = temp("data/fastq/from_cram/{SAMPLE}.2.fastq.gz"),
        fq0 = temp("data/fastq/from_cram/{SAMPLE}.0.fastq.gz")
    threads: 16
    resources:
        mem_mb = 32_000,
        runtime = 720
    params:
        tmp_prefix = lambda wildcards: f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/fastq/tmp/{wildcards.SAMPLE}",
    script:
       "bash_scripts/collate_fq.sh"
