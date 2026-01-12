import json
import pandas as pd

wildcard_constraints:
    ASSEMBLY = "GRCh38",
    CHROM = r"chr[0-9]{1,2}|chrX|chrY",
    SAMPLE = r"[0-9]{4,7}"


# map sample names to CRAM files
SMP2CRAM = {}
with open("json/cram_mapping.realigned.json") as f:
    dicts = json.load(f)
    for d in dicts:
        sample = d["ugrp_sample_id"]
        fh = d["cram_fh"]
        SMP2CRAM[sample] = fh


def get_cram_fh(wildcards):
    return SMP2CRAM[wildcards.SAMPLE]


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


rule call_snvs:
    input:
        ref = ELIFE_REF_FH,
        cram = get_cram_fh,
        cram_idx = lambda wildcards: get_cram_fh(wildcards) + ".crai",
        sif = "clara-parabricks_4.6.0-1.sif"
    output:
        vcf = "data/vcf/per-sample/{SAMPLE}.{ASSEMBLY}.vcf.gz",
        gvcf = "data/vcf/per-sample/{SAMPLE}.{ASSEMBLY}.g.vcf.gz",
    params:
        haploid_contigs_arg = lambda wildcards: "--haploid-contigs chrX,chrY" if SMP2SEX[wildcards.SAMPLE] == 1 else ""
    resources:
        mem_mb = 64_000,
        runtime = 60,
        gpu = 1,
        cpus_per_gpu = 64,
        slurm_account = "quinlan-gpu-rw",
        slurm_partition = "quinlan-gpu-rw",
        slurm_extra = "--exclusive"
    script:
        "bash_scripts/run_deepvariant_parabricks.sh"


# rule combine_trio_chrom_vcfs:
#     input:
#         sif = "glnexus_v1.2.7.sif",
#         gvcfs = expand("data/vcf/per-chrom/{SAMPLE}.{{ASSEMBLY}}.{{CHROM}}.g.vcf.gz", SAMPLE=list(set(all_samples))),
#     output: "data/vcf/per-chrom/{ASSEMBLY}.{CHROM}.joint_genotyped.vcf.gz"
#     threads: 8
#     params:
#         gl_nexus_prefix = lambda wildcards: f"gl_nexus_dbs/{wildcards.ASSEMBLY}_{wildcards.CHROM}"
#     resources:
#         mem_mb = 64_000
#     script:
#         "rules/bash_scripts/joint_genotype.sh"


# rule merge_trio_vcfs:
#     input:
#         vcfs = expand("data/vcf/per-chrom/{{ASSEMBLY}}.{CHROM}.joint_genotyped.vcf.gz", CHROM=CHROMS)
#     output: "data/vcf/merged/{ASSEMBLY}.joint_genotyped.vcf.gz"
#     threads: 16
#     resources:
#         runtime = 1440
#     shell:
#         """
#         module load bcftools
        
#         bcftools concat {input.vcfs} | bcftools view --threads {threads} | bgzip > {output}
#         """


# rule index_merged_vcf:
#     input: vcf = "data/vcf/merged/{ASSEMBLY}.joint_genotyped.vcf.gz"
#     output: "data/vcf/merged/{ASSEMBLY}.joint_genotyped.vcf.gz.tbi"
#     script:
#         "bash_scripts/index_vcf.sh"