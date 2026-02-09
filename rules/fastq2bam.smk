rule index_ref:
    input:
        reference = "data/ref/human_g1k_v38_decoy_phix.fasta",
        bwa_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/bwa-mem2-2.2.1_x64-linux/bwa-mem2"
    output:
        expand("data/ref/human_g1k_v38_decoy_phix.fasta.{suff}", suff=["pac", "ann", "amb"])
    shell:
        """
        {input.bwa_binary} -p data/ref/human_g1k_v38_decoy_phix {input.reference}
        """


rule fastq2bam:
    input:
        reference = "data/ref/human_g1k_v38_decoy_phix.fasta",
        fq1 = "data/fastq/{SAMPLE}.1.fastq.gz",
        fq2 = "data/fastq/{SAMPLE}.2.fastq.gz",
        bwa_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/bwa-mem2-2.2.1_x64-linux/bwa-mem2",
        sambamba_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/sambamba-1.0.1-linux-amd64-static"
    output:
        bam = temp("data/bam/{SAMPLE}.bam")
    threads: 32
    resources:
        mem_mb = 64_000
    shell:
        """        
        {input.bwa_binary} mem {input.reference} \
                {input.fq1} \
                {input.fq2} \
                -t {threads} | \
                {input.sambamba_binary} view -S -f bam -o {output.bam} /dev/stdin
        """


# NOTE: we can potentially use parabricks for alignment for a massive speedup.
# NOTE: parabricks will also markdups and sort in one go.
# rule fastq2bam_parabricks:
#     input:
#         reference = ELIFE_REF_FH,
#         fq1 = "data/fastq/{SAMPLE}.1.fastq.gz",
#         fq2 = "data/fastq/{SAMPLE}.2.fastq.gz",
#         sif = "clara-parabricks_4.6.0-1.sif",
#     output:
#         bam = temp("data/bam/{SAMPLE}.rg.sorted.bam"),
#     resources:
#         mem_mb = 64_000,
#         runtime = 60,
#         gpu = 1,
#         cpus_per_gpu = 64,
#         slurm_account = "quinlan-gpu-rw",
#         slurm_partition = "quinlan-gpu-rw",
#         slurm_extra = "--exclusive"
#     script:
#         "bash_scripts/run_fq2bam_parabricks.sh"


rule add_read_groups:
    input:
        bam = "data/bam/{SAMPLE}.bam"
    output:
        bam = temp("data/bam/{SAMPLE}.rg.bam")
    resources:
        mem_mb = 48_000
    shell:
        """
        module load gatk/4.6
        
        gatk --java-options "-Xmx48g" \
             AddOrReplaceReadGroups \
             -I {input.bam} \
             -O {output.bam} \
             --RGSM {wildcards.SAMPLE} \
             --RGLB lib1 \
             --RGPL ILLUMINA \
             --RGDS {wildcards.SAMPLE} \
             --RGPU unit1
        """


rule sort_bam:
    input:
        bam = "data/bam/{SAMPLE}.rg.bam",
        sambamba_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/sambamba-1.0.1-linux-amd64-static"
    output:
        bam = temp("data/bam/{SAMPLE}.rg.sorted.bam")
    threads: 8
    resources:
        mem_mb = 64_000
    params:
        tmpdir = "/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp/"
    shell:
        """       
        {input.sambamba_binary} sort \
                      -t {threads} \
                      -m 32G \
                      -p \
                      -o {output.bam} \
                      --tmpdir {params.tmpdir} \
                      {input.bam}
        """


rule index_bam:
    input:
        bam = "data/bam/{SAMPLE}.rg.sorted.bam",
        sambamba_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/sambamba-1.0.1-linux-amd64-static"
    output:
        idx = "data/bam/{SAMPLE}.rg.sorted.bam.bai"
    threads: 8
    shell:
        """        
        {input.sambamba_binary} index -t {threads} {input.bam}
        """ 

    
rule convert_to_cram:
    input:
        bam = "data/bam/{SAMPLE}.rg.sorted.bam",
        reference = ELIFE_REF_FH,
    output:
        cram = temp("data/cram/{SAMPLE}.cram")
    threads: 8
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

# NOTE: careful, this removes the embedded reference!
rule mark_duplicates:
    input:
        cram = "data/cram/{SAMPLE}.cram",
        reference = "data/ref/human_g1k_v38_decoy_phix.fasta",
    output:
        cram = "data/cram/{SAMPLE}.dupmarked.cram",
        metrics = "data/markdup_metrics/{SAMPLE}.metrics.txt"
    params:
        tmpdir = "/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp/"
    resources:
        mem_mb = 64_000,
        runtime = 1_080
    shell:
        """
        module load gatk/4.6

        gatk --java-options "-Xmx48g" \
             MarkDuplicates \
             -I {input.cram} \
             -O {output.cram} \
             -R {input.reference} \
             --REMOVE_DUPLICATES false \
             --TMP_DIR {params.tmpdir} \
             --METRICS_FILE {output.metrics} \
             --VALIDATION_STRINGENCY SILENT
        """


rule index_cram:
    input:
        cram = "data/cram/{SAMPLE}.dupmarked.cram",
    output:
        idx = "data/cram/{SAMPLE}.dupmarked.cram.crai"
    threads: 8
    resources:
        mem_mb = 32_000
    shell:
        """
        module load samtools
        
        samtools index -@ {threads} {input.cram}
        """

