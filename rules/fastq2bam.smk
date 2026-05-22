rule fastq2bam:
    input:
        ref = config["alignment_ref"],
        fq1 = lambda wildcards: f"data/fastq/{SAMPLE2SOURCE[wildcards.SAMPLE]}/{wildcards.SAMPLE}.1.fastq.gz",
        fq2 = lambda wildcards: f"data/fastq/{SAMPLE2SOURCE[wildcards.SAMPLE]}/{wildcards.SAMPLE}.2.fastq.gz",
    output:
        bam = temp("data/bam/{SAMPLE}.bam")
    threads: 32
    resources:
        runtime = 360,
        mem_mb = 64_000
    script:
        "bash_scripts.align.sh"


# NOTE: we can potentially use parabricks for alignment for a massive speedup.
# parabricks will also markdups and sort in one go.
# rule fastq2bam_parabricks:
#     input:
#         reference = config["alignment_ref"],
#         fq1 = lambda wildcards: f"data/fastq/{SAMPLE2SOURCE[wildcards.SAMPLE]}/{wildcards.SAMPLE}.1.fastq.gz",
#         fq2 = lambda wildcards: f"data/fastq/{SAMPLE2SOURCE[wildcards.SAMPLE]}/{wildcards.SAMPLE}.2.fastq.gz",
#         sif = "clara-parabricks_4.6.0-1.sif",
#     output:
#         bam = temp("data/bam/{SAMPLE}.dupmarked.bam"),
#         metrics = "data/markdup_metrics/{SAMPLE}.metrics.txt"
#     params:
#         memory_limit = 32_000
#     resources:
#         mem_mb = 64_000,
#         runtime = 90,
#         gpu = 1,
#         cpus_per_gpu = 64,
#         slurm_account = "quinlan-gpu-rw",
#         slurm_partition = "quinlan-gpu-rw",
#         slurm_extra = "--exclusive"
#     script:
#         "bash_scripts/run_fq2bam_parabricks.sh"


# rule add_read_groups:
#     input:
#         bam = "data/bam/{SAMPLE}.bam"
#     output:
#         bam = temp("data/bam/{SAMPLE}.rg.bam")
#     resources:
#         mem_mb = 48_000
#     shell:
#         """
#         module load gatk/4.6
        
#         gatk --java-options "-Xmx48g" \
#              AddOrReplaceReadGroups \
#              -I {input.bam} \
#              -O {output.bam} \
#              --RGSM {wildcards.SAMPLE} \
#              --RGLB lib1 \
#              --RGPL ILLUMINA \
#              --RGDS {wildcards.SAMPLE} \
#              --RGPU unit1
#         """


# rule sort_bam:
#     input:
#         bam = "data/bam/{SAMPLE}.rg.bam"
#     output:
#         bam = temp("data/bam/{SAMPLE}.rg.sorted.bam")
#     threads: 8
#     resources:
#         mem_mb = 64_000
#     params:
#         tmpdir = "/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp/"
#     shell:
#         """   
#         module load sambamba/040218    
#         sambamba sort \
#                       -t {threads} \
#                       -m 32G \
#                       -p \
#                       -o {output.bam} \
#                       --tmpdir {params.tmpdir} \
#                       {input.bam}
#         """


# rule index_bam:
#     input:
#         bam = "data/bam/{SAMPLE}.rg.sorted.bam",
#         sambamba_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/sambamba-1.0.1-linux-amd64-static"
#     output:
#         idx = "data/bam/{SAMPLE}.rg.sorted.bam.bai"
#     threads: 8
#     shell:
#         """        
#         {input.sambamba_binary} index -t {threads} {input.bam}
#         """ 


rule mark_duplicates:
    input:
        bam = "data/bam/{SAMPLE}.bam",
        reference = config["alignment_ref"],
    output:
        bam = "data/bam/{SAMPLE}.dupmarked.bam",
        metrics = "data/markdup_metrics/{SAMPLE}.metrics.txt"
    params:
        tmpdir = "/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp/"
    resources:
        mem_mb = 48_000,
        runtime = 360
    shell:
        """
        module load gatk/4.6

        gatk --java-options "-Xmx32g" \
             MarkDuplicates \
             -I {input.bam} \
             -O {output.bam} \
             -R {input.reference} \
             --REMOVE_DUPLICATES false \
             --TMP_DIR {params.tmpdir} \
             --METRICS_FILE {output.metrics} \
             --VALIDATION_STRINGENCY SILENT
        """


rule convert_to_cram:
    input:
        bam = "data/bam/{SAMPLE}.dupmarked.bam",
        reference = config["bam_to_cram_ref"],
    output:
        cram = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/cram/{SAMPLE}.cram"
    threads: 8
    resources:
        runtime = 180
    shell:
        """
        module load samtools

        samtools view -@ {threads} \
                      -C \
                      -o {output.cram} \
                      -T {input.reference} \
                      --output-fmt-option embed_ref=1 \
                      {input.bam}
        """



rule index_cram:
    input:
        cram = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/cram/{SAMPLE}.cram",
    output:
        idx = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/cram/{SAMPLE}.cram.crai"
    threads: 8
    resources:
        mem_mb = 32_000
    shell:
        """
        module load samtools
        
        samtools index -@ {threads} {input.cram}
        """

