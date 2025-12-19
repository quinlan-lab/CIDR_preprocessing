ELIFE_REF_FH = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"


rule fastq2bam:
    input:
        reference = ELIFE_REF_FH,
        fq1 = "data/fastq/{SAMPLE}.1.clean.fastq.gz",
        fq2 = "data/fastq/{SAMPLE}.2.clean.fastq.gz"
    output:
        bam = temp("data/bam/{SAMPLE}.bam")
    threads: 16
    resources:
        mem_mb = 32_000
    shell:
        """
        module load bwa samtools
        
        bwa mem {input.reference} \
                {input.fq1} \
                {input.fq2} \
                -t {threads} | \
                samtools sort -o {output.bam}
        """


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
        
        gatk -Xmx48g \
             AddOrReplaceReadGroups \
             -I {input.bam} \
             -O {output.bam} \
             --RGSM {wildcards.SAMPLE} \
             --RGLB lib1 \
             --RGPL illumina \
             --RGDS {wildcards.SAMPLE} \
             --RGPU unit1
        """


rule sort_bam:
    input:
        bam = "data/bam/{SAMPLE}.rg.bam"
    output:
        bam = temp("data/bam/{SAMPLE}.rg.sorted.bam")
    threads: 8
    resources:
        mem_mb = 16_000
    params:
        tmpdir = "/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp/"
    shell:
        """
        module load samtools
        
        samtools sort -@ {threads} \
                      -Ob \
                      -o {output.bam} \
                      -T {params.tmpdir} \
                      -m 16G \
                      {input.bam}
        """


rule index_bam:
    input:
        bam = "data/bam/{SAMPLE}.rg.sorted.bam"
    output:
        idx = "data/bam/{SAMPLE}.rg.sorted.bam.bai"
    threads: 8
    shell:
        """
        module load samtools
        
        samtools index -@ {threads} {input.bam}
        """ 

    
rule convert_to_cram:
    input:
        bam = "data/bam/{SAMPLE}.rg.sorted.bam",
        reference = ELIFE_REF_FH,
    output:
        cram = "data/cram/{SAMPLE}.cram"
    shell:
        """
        module load samtools

        samtools view -@ {threads} \
                      -C \
                      -o {output.cram} \
                      -T {input.reference} \
                      --output-fmt-option embed_ref=1
        """


rule mark_duplicates:
    input:
        cram = "data/cram/{SAMPLE}.cram"
    output:
        cram = "data/cram/{SAMPLE}.dupmarked.cram",
        idx = "data/cram/{SAMPLE}.dupmarked.cram.crai"
    params:
        tmpdir = "/scratch/ucgd/lustre-labs/quinlan/u1006375/samtools_tmp/"
    resources:
        mem_mb = 48_000
    shell:
        """
        module load gatk/4.6

        gatk -Xmx48g \
             MarkDuplicates \
             -I {input.cram} \
             -O {output.cram} \
             --REMOVE_DUPLICATES false \
             --CREATE_INDEX true \
             --TMP_DIR {params.tmpdir} \
             --VALIDATION_STRINGENCY SILENT
        """

