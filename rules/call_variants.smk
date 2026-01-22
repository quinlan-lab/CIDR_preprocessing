import json
import pandas as pd
import numpy as np

wildcard_constraints:
    ASSEMBLY = "GRCh38",
    CHROM = r"chr[0-9]{1,2}|chrX|chrY",
    SAMPLE = r"[0-9]{4,7}"


quinlan_chroms = [f"chr{c}" for c in range(1, 12)]
ucgd_chroms = [f"chr{c}" for c in range(12, 23)] + ["chrX", "chrY"]


gpu_chroms = [f"chr{c}" for c in range(1, 12)]
cpu_chroms = [f"chr{c}" for c in range(12, 23)] + ["chrX", "chrY"]

quinlan_chroms = cpu_chroms[:6]
ucgd_chroms = cpu_chroms[6:]

# map sample names to CRAM files
smp2cram = {}
with open("/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/json/cram_mapping.realigned.json") as f:
    dicts = json.load(f)
    for d in dicts:
        sample = d["ugrp_sample_id"]
        fh = d["cram_fh"]
        smp2cram[sample] = fh

def get_cram_fh(wildcards):
    return smp2cram[wildcards.SAMPLE]
def get_cram_idx_fh(wildcards):
    return smp2cram[wildcards.SAMPLE] + ".crai"

ELIFE_REF_FH = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"

MASTER_PED = "/scratch/ucgd/lustre-labs/quinlan/u0890814/CIDR_4Gen/ped_files/master_ped_all_info.ped"
# get sample IDs in new CIDR CEPH cohort
sample_info = pd.read_csv(
    MASTER_PED,
    sep="\t",
    dtype={"UGRP_Lab_ID": str, "Gender": str}
)
sample_info = sample_info[sample_info["Gender"].isin(["1", "2"])]
SMP2SEX = dict(zip(sample_info["UGRP_Lab_ID"], list(map(int, sample_info["Gender"]))))


# rule call_snvs_parabricks:
#     input:
#         ref = ELIFE_REF_FH,
#         cram = get_cram_fh,
#         cram_idx = get_cram_idx_fh,
#         sif = "clara-parabricks_4.6.0-1.sif"
#     output:
#         vcf = "data/vcf/per-sample/{SAMPLE}.{ASSEMBLY}.vcf.gz",
#         gvcf = "data/vcf/per-sample/{SAMPLE}.{ASSEMBLY}.g.vcf.gz",
#     params:
#         haploid_contigs_arg = lambda wildcards: "--haploid-contigs chrX,chrY" if SMP2SEX[wildcards.SAMPLE] == 1 else ""
#     resources:
#         mem_mb = 64_000,
#         runtime = 60,
#         gpu = 1,
#         cpus_per_gpu = 64,
#         slurm_account = "quinlan-gpu-rw",
#         slurm_partition = "quinlan-gpu-rw",
#         slurm_extra = "--exclusive"
#     script:
#         "bash_scripts/run_deepvariant_parabricks.sh"


# rule call_snvs_cpu:
#     input:
#         ref = ELIFE_REF_FH,
#         cram = get_cram_fh,
#         cram_idx = get_cram_idx_fh,
#         sif = "deepvariant_1.9.0.sif"
#     output:
#         vcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf.gz",
#         gvcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.g.vcf.gz",
#     params:
#         haploid_contigs_arg = lambda wildcards: "--haploid_contigs chrX,chrY" if SMP2SEX[wildcards.SAMPLE] == 1 else ""
#     resources:
#         mem_mb = 32_000,
#         runtime = 360,
#         slurm_account = lambda wildcards: "ucgd-rw" if wildcards.CHROM in ucgd_chroms else "quinlan-rw",
#         slurm_partition = lambda wildcards: "ucgd-rw" if wildcards.CHROM in ucgd_chroms else "quinlan-rw"
#     threads: 16
#     script:
#         "bash_scripts/run_deepvariant.sh"



def get_haploid_contigs_arg(wildcards, cpu_chroms):
    if SMP2SEX[wildcards.SAMPLE] == 2:
        return ""
    else:
        if wildcards.CHROM in cpu_chroms:
            return "--haploid_contigs chrX,chrY"
        else:
            return "--haploid-contigs chrX,chrY"

rule call_snvs_cpu_gpu:
    input:
        ref = ELIFE_REF_FH,
        cram = get_cram_fh,
        cram_idx = get_cram_idx_fh,
    output:
        vcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf.gz",
        gvcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.g.vcf.gz",
    params:
        sif = lambda wildcards: "deepvariant_1.9.0.sif" if wildcards.CHROM in cpu_chroms else "clara-parabricks_4.6.0-1.sif",
        haploid_contigs_arg = lambda wildcards: get_haploid_contigs_arg(wildcards, cpu_chroms),
        gpu = lambda wildcards: 0 if wildcards.CHROM in cpu_chroms else 1
    resources:
        mem_mb = lambda wildcards: 64_000 if wildcards.CHROM in gpu_chroms else 32_000,
        runtime = lambda wildcards: 60 if wildcards.CHROM in gpu_chroms else 360,
        gpu = lambda wildcards: 1 if wildcards.CHROM in gpu_chroms else None,
        cpus_per_gpu = lambda wildcards: 64 if wildcards.CHROM in gpu_chroms else None,
        slurm_account = lambda wildcards: "ucgd-rw" if wildcards.CHROM in ucgd_chroms else "quinlan-rw" if wildcards.CHROM in quinlan_chroms else "quinlan-gpu-rw",
        slurm_partition = lambda wildcards: "ucgd-rw" if wildcards.CHROM in ucgd_chroms else "quinlan-rw" if wildcards.CHROM in quinlan_chroms else "quinlan-gpu-rw",
        slurm_extra = lambda wildcards: "--exclusive" if wildcards.CHROM in gpu_chroms else None
    threads: 16
    script:
        "bash_scripts/run_deepvariant_cpu_gpu.sh"

# rule merge_chrom_vcfs:
#     input:
#         gvcfs = expand("data/vcf/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{CHROM}.g.vcf.gz", CHROM=CHROMS)
#     output: "data/vcf/per-sample/{SAMPLE}.{ASSEMBLY}.vcf.gz"
#     threads: 16
#     resources:
#         runtime = 1440
#     shell:
#         """
#         module load bcftools
        
#         bcftools concat {input.gvcfs} | bcftools view --threads {threads} | bgzip > {output}
#         """


