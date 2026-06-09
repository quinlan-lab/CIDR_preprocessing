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
        runtime = 240,
        slurm_account = "ucgd-rw",#lambda wildcards: accts[int(wildcards.SAMPLE) % 2],
        slurm_partition = "ucgd-rw",#lambda wildcards: accts[int(wildcards.SAMPLE) % 2]
    threads: 24
    script:
        "bash_scripts/run_deepvariant.sh"


# def get_haploid_contigs_arg(wildcards, cpu_chroms):
#     if SMP2SEX[wildcards.SAMPLE] == 2:
#         return ""
#     else:
#         if wildcards.CHROM in cpu_chroms:
#             return "--haploid_contigs chrX,chrY"
#         else:
#             return "--haploid-contigs chrX,chrY"


# # NOTE: this rule supersedes above two rules -- bash script contains
# # logic for running DV on either CPU or GPU, depending on the chrom
# # we're processing.
# rule call_snvs_cpu_gpu:
#     input:
#         ref = ELIFE_REF_FH,
#         cram = get_new_cram_fh,
#         cram_idx = get_new_cram_idx_fh,
#     output:
#         vcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf.gz",
#         gvcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.g.vcf.gz",
#     params:
#         sif = lambda wildcards: "deepvariant_1.9.0.sif" if wildcards.CHROM in cpu_chroms else "clara-parabricks_4.6.0-1.sif",
#         haploid_contigs_arg = lambda wildcards: get_haploid_contigs_arg(wildcards, cpu_chroms),
#         gpu = lambda wildcards: 0 if wildcards.CHROM in cpu_chroms else 1
#     resources:
#         mem_mb = lambda wildcards: 64_000 if wildcards.CHROM in gpu_chroms else 32_000,
#         runtime = lambda wildcards: 60 if wildcards.CHROM in gpu_chroms else 360,
#         gpu = lambda wildcards: 1 if wildcards.CHROM in gpu_chroms else None,
#         cpus_per_gpu = lambda wildcards: 64 if wildcards.CHROM in gpu_chroms else None,
#         slurm_account = lambda wildcards: "ucgd-rw" if wildcards.CHROM in ucgd_chroms else "quinlan-rw" if wildcards.CHROM in quinlan_chroms else "quinlan-gpu-rw",
#         slurm_partition = lambda wildcards: "ucgd-rw" if wildcards.CHROM in ucgd_chroms else "quinlan-rw" if wildcards.CHROM in quinlan_chroms else "quinlan-gpu-rw",
#         slurm_extra = lambda wildcards: "--exclusive" if wildcards.CHROM in gpu_chroms else "--exclude=rw159"
#     threads: 16
#     script:
#         "bash_scripts/run_deepvariant_cpu_gpu.sh"


