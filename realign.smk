import json
import pandas as pd

CHROMS = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]

# we use the DRAGEN REF as the basis for extracting FASTQ for CIDR samples
DRAGEN_REF_FH = "/scratch/ucgd/lustre-labs/quinlan/u6070793/master_files/hg38.fa"
# we re-align CIDR to the HG38 version used for the ELIFE samples
ELIFE_REF_FH = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"

PREF = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH"

PROVENANCE = pd.read_csv("PROVENANCE.tsv", sep="\t", dtype={"UGRP_Lab_ID": str})

MASTER_PED = "/scratch/ucgd/lustre-labs/quinlan/u0890814/CIDR_4Gen/ped_files/master_ped_all_info.ped"
# get sample info for full cohort from Julia's "ped"
SAMPLE_INFO = pd.read_csv(
    MASTER_PED,
    sep="\t",
    dtype={"UGRP_Lab_ID": str, "Gender": str}
).dropna(subset=["Sequencing"])
SAMPLE_INFO = SAMPLE_INFO[SAMPLE_INFO["Gender"].isin(["1", "2"])]
SMP2SEX = dict(zip(SAMPLE_INFO["UGRP_Lab_ID"], list(map(int, SAMPLE_INFO["Gender"]))))


def get_new_cram_fh(wildcards):
    sample_map =  FILE_MAP[FILE_MAP["UGRP_Lab_ID"] == wildcards.SAMPLE]
    assert sample_map.shape[0] == 1
    return sample_map["final_cram_path"].to_list().pop()

def get_orig_cram_fh(wildcards):
    sample_map =  FILE_MAP[FILE_MAP["UGRP_Lab_ID"] == wildcards.SAMPLE]
    assert sample_map.shape[0] == 1
    return sample_map["orig_cram_path"].to_list().pop()

def get_new_cram_idx_fh(wildcards):
    return get_new_cram_fh(wildcards) + ".crai"
    

include: "rules/call_variants.smk"
include: "rules/cram2fastq.smk"
include: "rules/fastq2bam.smk"


wildcard_constraints:
    ASSEMBLY = "GRCh38",
    CHROM = r"chr[0-9]{1,2}|chrX|chrY",
    SAMPLE = r"[0-9]{4,7}"


samples = PROVENANCE[PROVENANCE["sequencing_provenance"] == "UofU_rd2,eLife"]["UGRP_Lab_ID"].to_list()

rule all:
    input:
        # expand("data/vcf/per-chrom/{SAMPLE}.GRCh38.{CHROM}.g.vcf.gz", SAMPLE = ["300100"], CHROM=["chr21"]),
        expand("data/cram/{SAMPLE}.dupmarked.cram.crai", SAMPLE = samples)


# rule joint_genotype:
#     input:
#         sif = "glnexus_v1.2.7.sif",
#         gvcfs = expand("data/vcf/per-chrom/{SAMPLE}.{{ASSEMBLY}}.{{CHROM}}.g.vcf.gz", SAMPLE=all_samples),
#     output: "data/vcf/joint_genotyped/{ASSEMBLY}.{CHROM}.joint_genotyped.bcf"
#     threads: 16
#     params:
#         gl_nexus_prefix = lambda wildcards: f"gl_nexus_dbs/{wildcards.ASSEMBLY}_{wildcards.CHROM}"
#     resources:
#         mem_mb = 64_000
#     script:
#         "rules/bash_scripts/joint_genotype.sh"