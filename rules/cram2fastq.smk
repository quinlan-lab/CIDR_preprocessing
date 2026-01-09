import json


# we use the DRAGEN REF as the basis for extracting FASTQ for CIDR samples
DRAGEN_REF_FH = "/scratch/ucgd/lustre-labs/quinlan/u6070793/master_files/hg38.fa"


# map sample names to CRAM files
SMP2CRAM = {}
with open("json/cram_mapping.json") as f:
    dicts = json.load(f)
    for d in dicts:
        sample = d["sample"]
        fh = d["cram_fh"]
        SMP2CRAM[sample] = fh


def get_cram_fh(wildcards):
    return SMP2CRAM[wildcards.SAMPLE]


rule namesort_cram:
    input:
        cram = get_cram_fh,
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
        reference = DRAGEN_REF_FH,
    output:
        fq1 = temp("data/fastq/{SAMPLE}.1.fastq.gz"),
        fq2 = temp("data/fastq/{SAMPLE}.2.fastq.gz")
    threads: 8
    resources:
        mem_mb = 32_000
    shell:
        """
        module load samtools
        
        samtools fastq -@ {threads} \
                       --reference {input.reference} \
                       -1 {output.fq1} \
                       -2 {output.fq2} \
                       {input.cram}
        """


rule clean_with_fastq:
    input:
        fq1 = "data/fastq/{SAMPLE}.1.fastq.gz",
        fq2 = "data/fastq/{SAMPLE}.2.fastq.gz",
    output:
        fq1_clean = temp("data/fastq/{SAMPLE}.1.clean.fastq.gz"),
        fq2_clean = temp("data/fastq/{SAMPLE}.2.clean.fastq.gz"),
        html_report = "reports/{SAMPLE}-fastp-report.html",
        json_report = "reports/{SAMPLE}-fastp-report.json"
    threads: 16
    log: "logs/fastp/{SAMPLE}.log"
    shell:
        """
        module load fastp

        fastp --in1 {input.fq1} \
                --in2 {input.fq2} \
                --out1 {output.fq1_clean} \
                --out2 {output.fq2_clean} \
                --thread {threads} \
                --disable_quality_filtering \
                --disable_adapter_trimming \
                --disable_trim_poly_g \
                --html {output.html_report} \
                --json {output.json_report}
        """
