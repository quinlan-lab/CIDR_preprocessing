import json

DRAGEN_REF_FH = "/scratch/ucgd/lustre-labs/quinlan/u6070793/hg38.fa"


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
        cram = get_cram_fh
    output:
        cram = temp("data/cram/{SAMPLE}.sorted.cram")
    shell:
        """
        module load samtools
        
        samtools sort -@ {threads} \
                      -n \
                      -O CRAM \
                      -o {output.cram} \
                      {input.cram}
        """


rule cram2fastq:
    input:
        cram = get_cram_fh,
        reference = DRAGEN_REF_FH,
    output:
        fq1 = temp("data/fastq/{SAMPLE}.1.fastq.gz"),
        fq2 = temp("data/fastq/{SAMPLE}.2.fastq.gz")
    threads: 8
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
        fq1_clean = temp("data/fastq/{SAMPLE}.1.clean.fastq"),
        fq2_clean = temp("data/fastq/{SAMPLE}.2.clean.fastq"),
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


rule compress_fastq:
    input:
        fq = "data/fastq/{SAMPLE}.{PAIR}.clean.fastq"
    output:
        fq = "data/fastq/{SAMPLE}.{PAIR}.clean.fastq.gz"
    shell:
        """
        module load bgzip
        
        bgzip {input.fq}
        """
