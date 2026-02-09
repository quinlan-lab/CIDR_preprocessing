import json
import pandas as pd

CHROMS = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]
gpu_chroms = CHROMS[:2]
cpu_chroms = CHROMS[2:]
quinlan_chroms = cpu_chroms[10:]
ucgd_chroms = cpu_chroms[:10]

# we use the DRAGEN REF as the basis for extracting FASTQ for CIDR samples
DRAGEN_REF_FH = "/scratch/ucgd/lustre-labs/quinlan/u6070793/master_files/hg38.fa"
# we re-align CIDR to the HG38 version used for the ELIFE samples
ELIFE_REF_FH = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"

# map sample names to original CRAM files, from which we will
# extract FASTQ if necessary
SMP2CRAM_ORIG = {}
with open("json/cram_mapping.json") as f:
    dicts = json.load(f)
    for d in dicts:
        sample = d["sample"]
        fh = d["cram_fh"]
        SMP2CRAM_ORIG[sample] = fh

def get_orig_cram_fh(wildcards):
    return SMP2CRAM_ORIG[wildcards.SAMPLE]

# map sample names to new CRAM file handles after realignment, which 
# we'll use to call variants
SMP2CRAM_NEW = {}
with open("/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/json/cram_mapping.realigned.json") as f:
    dicts = json.load(f)
    for d in dicts:
        sample = d["ugrp_sample_id"]
        fh = d["cram_fh"]
        SMP2CRAM_NEW[sample] = fh


def get_new_cram_fh(wildcards):
    return SMP2CRAM_NEW[wildcards.SAMPLE]

def get_new_cram_idx_fh(wildcards):
    return SMP2CRAM_NEW[wildcards.SAMPLE] + ".crai"


MASTER_PED = "/scratch/ucgd/lustre-labs/quinlan/u0890814/CIDR_4Gen/ped_files/master_ped_all_info.ped"
# get sample IDs in new CIDR CEPH cohort
sample_info = pd.read_csv(
    MASTER_PED,
    sep="\t",
    dtype={"UGRP_Lab_ID": str, "Gender": str}
)
sample_info = sample_info[sample_info["Gender"].isin(["1", "2"])]
SMP2SEX = dict(zip(sample_info["UGRP_Lab_ID"], list(map(int, sample_info["Gender"]))))


include: "rules/call_variants.smk"
include: "rules/cram2fastq.smk"
include: "rules/fastq2bam.smk"


wildcard_constraints:
    ASSEMBLY = "GRCh38",
    CHROM = r"chr[0-9]{1,2}|chrX|chrY",
    SAMPLE = r"[0-9]{4,7}"


MASTER_PED = "/scratch/ucgd/lustre-labs/quinlan/u0890814/CIDR_4Gen/ped_files/master_ped_all_info.ped"

# get sample IDs in new CIDR CEPH cohort
sample_info = pd.read_csv(
    MASTER_PED,
    sep="\t",
    dtype={"UGRP_Lab_ID": str}
)

cidr = sample_info[sample_info["Sequencing"].isin(["CIDR-Illumina_short-read", "CIDR-Illumina_short-read_top-up"])]["UGRP_Lab_ID"].to_list()
elife = sample_info[sample_info["Sequencing"] == "WashU-Illumina_short-read"]["UGRP_Lab_ID"].to_list()
scott = sample_info[sample_info["Sequencing"] == "2025UofU"]["UGRP_Lab_ID"].to_list()

scott = scott + ["300010", "200093", "200112", "200138"]
scott = [s for s in scott if s not in ("210071", "130053", "300010")] # no sex for 210071


rule all:
    input:
        # expand("data/vcf/joint_genotyped/{ASSEMBLY}.{CHROM}.joint_genotyped.bcf", ASSEMBLY=["GRCh38"], CHROM=CHROMS)
        expand("data/vcf/per-chrom/{SAMPLE}.GRCh38.{CHROM}.g.vcf.gz", SAMPLE = scott, CHROM=CHROMS),
        # expand("data/cram/{SAMPLE}.dupmarked.cram.crai", SAMPLE = scott)


rule joint_genotype:
    input:
        sif = "glnexus_v1.2.7.sif",
        gvcfs = expand("data/vcf/per-chrom/{SAMPLE}.{{ASSEMBLY}}.{{CHROM}}.g.vcf.gz", SAMPLE=list(set(cidr + elife))),
    output: "data/vcf/joint_genotyped/{ASSEMBLY}.{CHROM}.joint_genotyped.bcf"
    threads: 16
    params:
        gl_nexus_prefix = lambda wildcards: f"gl_nexus_dbs/{wildcards.ASSEMBLY}_{wildcards.CHROM}"
    resources:
        mem_mb = 64_000
    script:
        "rules/bash_scripts/joint_genotype.sh"