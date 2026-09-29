rule decompress:
    input:
        orad_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/orad.2.7.0.linux/orad",
        ora_list = lambda wildcards: [str(s) for s in sample2ora[wildcards.SAMPLE]],
    output:
        placeholder_fastq = "data/fastq/{SAMPLE}_0.fq.gz_1",
    params:
        n_input_files = lambda wildcards: len(sample2ora[wildcards.SAMPLE]) - 1,
        prefix = PREF
    threads: 16
    script:
        "bash_scripts/decompress_ora.sh"


rule concatenate:
    input:
        placeholder_fastq = "data/fastq/{SAMPLE}_0.fq.gz_1"
    output:
        fq1 = temp("data/fastq/from_ora/{SAMPLE}.1.fastq.gz"),
        fq2 = temp("data/fastq/from_ora/{SAMPLE}.2.fastq.gz")
    params:
        n_input_files = lambda wildcards: len(sample2ora[wildcards.SAMPLE]),
        fastq1_list = lambda wildcards: [f"data/fastq/{wildcards.SAMPLE}_{n_fq}.fq.gz_1" for n_fq in range(len(sample2ora[wildcards.SAMPLE]))],
        fastq2_list = lambda wildcards: [f"data/fastq/{wildcards.SAMPLE}_{n_fq}.fq.gz_2" for n_fq in range(len(sample2ora[wildcards.SAMPLE]))]
    script:
        "bash_scripts/combine_fq.sh"