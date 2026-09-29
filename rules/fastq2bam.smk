ACCT_OPTIONS = ["quinlan-rw", "ucgd-rw"]

rule fastq2bam:
    input:
        alignment_ref = config["alignment_ref"],
        cram_ref = config["bam_to_cram_ref"],
        alt_contigs = config["alignment_alt_contigs"],
        samblaster_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/samblaster",
        postalt = "/uufs/chpc.utah.edu/sys/installdir/bwa/0.7.17/bwakit/bwa-postalt.js",
        fq1 = lambda wildcards: f"data/fastq/{SAMPLE2SOURCE[wildcards.SAMPLE]}/{wildcards.SAMPLE}.1.fastq.gz",
        fq2 = lambda wildcards: f"data/fastq/{SAMPLE2SOURCE[wildcards.SAMPLE]}/{wildcards.SAMPLE}.2.fastq.gz",
    output:
        cram = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/cram/{SAMPLE}.cram",
        cram_idx = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/cram/{SAMPLE}.cram.crai"
    threads: 32
    params:
        K = 96_000_000,
        read_group = lambda wildcards: f"@RG\\tID:{wildcards.SAMPLE}\\tLB:lib1\\tPL:ILLUMINA\\tSM:{wildcards.SAMPLE}\\tPU:{wildcards.SAMPLE}",
        tmp_prefix = lambda wildcards: f"/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp/{wildcards.SAMPLE}"
    resources:
        runtime = 960,
        mem_mb = 64_000,
        slurm_account = "ucgd-rw", #lambda wildcards: ACCT_OPTIONS[int(wildcards.SAMPLE) % 2],
        slurm_partition = "ucgd-rw", #lambda wildcards: ACCT_OPTIONS[int(wildcards.SAMPLE) % 2],
    script:
        "bash_scripts/align.sh"