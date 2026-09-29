ACCT_OPTIONS = ["quinlan-rw", "ucgd-rw"]

rule call_snvs_cpu:
    input:
        ref = config["bam_to_cram_ref"],
        cram = get_new_cram_fh,
        cram_idx = get_new_cram_idx_fh,
        sif = "deepvariant_1.9.0.sif"
    output:
        vcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf.gz",
        gvcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.g.vcf.gz",
    params:
        haploid_contigs_arg = lambda wildcards: "--haploid_contigs chrX,chrY" if SMP2SEX[wildcards.SAMPLE] == 1 else ""
    resources:
        mem_mb = 32_000,
        runtime = 960,
        slurm_account = lambda wildcards: ACCT_OPTIONS[int(wildcards.SAMPLE) % 2],
        slurm_partition = lambda wildcards: ACCT_OPTIONS[int(wildcards.SAMPLE) % 2],
    threads: 16
    script:
        "bash_scripts/run_deepvariant.sh"

GPU_OPTIONS = ["quinlan-gpu-rw", "ucgd-gpu-rw"]

rule call_snvs_gpu:
    input:
        ref = config["bam_to_cram_ref"],
        cram = get_new_cram_fh,
        cram_idx = get_new_cram_idx_fh,
        sif = "clara-parabricks_4.7.1-1.sif"
    output:
        vcf = "data/vcf/per-sample/{SAMPLE}.{ASSEMBLY}.vcf.gz",
        gvcf = "data/vcf/per-sample/{SAMPLE}.{ASSEMBLY}.g.vcf.gz",
    params:
        haploid_contigs_arg = lambda wildcards: "--haploid-contigs chrX,chrY" if SMP2SEX[wildcards.SAMPLE] == 1 else ""
    threads: 16
    resources:
        mem_mb = 64_000,
        runtime = 60,
        cpus_per_gpu = 64,
        gpu = 1,
        slurm_account = lambda wildcards: GPU_OPTIONS[int(wildcards.SAMPLE) % 2],
        slurm_partition = lambda wildcards: GPU_OPTIONS[int(wildcards.SAMPLE) % 2],
        slurm_extra = "--exclusive"
    script:
        "bash_scripts/run_deepvariant_parabricks.sh"


# rule split_by_chrom:
#     input:
#         vcf = "data/vcf/per-sample/{SAMPLE}.{ASSEMBLY}.vcf.gz",
#         gvcf = "data/vcf/per-sample/{SAMPLE}.{ASSEMBLY}.g.vcf.gz",
#     output:
#         vcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf.gz",
#         gvcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.g.vcf.gz",
#     threads: 4
#     shell:
#         """
#         module load bcftools
        
#         bcftools view --threads {threads} -Ob -o {output.vcf} {input.vcf} {wildcards.CHROM}
#         bcftools view --threads {threads} -Ob -o {output.gvcf} {input.gvcf} {wildcards.CHROM}
#         """
        

