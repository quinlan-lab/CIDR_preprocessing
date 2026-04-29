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

SMP2CRAM_ORIG = {}
with open("json/cram_mapping.to_extract.json") as f:
    dicts = json.load(f)
    for d in dicts:
        sample = d["ugrp_sample_id"]
        fh = d["cram_fh"]
        SMP2CRAM_ORIG[sample] = fh

def get_orig_cram_fh(wildcards):
    return SMP2CRAM_ORIG[wildcards.SAMPLE]

# map sample names to new CRAM file handles after realignment, which 
# we'll use to call variants
SMP2CRAM_NEW = {}
with open("json/cram_mapping.realigned.json") as f:
    dicts = json.load(f)
    for d in dicts:
        sample = d["ugrp_sample_id"]
        fh = d["cram_fh"]
        SMP2CRAM_NEW[sample] = fh


def get_new_cram_fh(wildcards):
    return SMP2CRAM_NEW[wildcards.SAMPLE]

def get_new_cram_idx_fh(wildcards):
    return SMP2CRAM_NEW[wildcards.SAMPLE] + ".crai"


# def get_fastq_fh(wildcards):


include: "rules/call_variants.smk"

### NOTE: comment out this rule if we're running on the "alignment" (rather)
# than "realignment" samples
include: "rules/cram2fastq.smk"
include: "rules/fastq2bam.smk"


wildcard_constraints:
    ASSEMBLY = "GRCh38",
    CHROM = r"chr[0-9]{1,2}|chrX|chrY",
    SAMPLE = r"[0-9]{4,7}"


### NOTE: these are samples for whom we do nothing!
elife_provenance = (["eLife"])
samples_to_keep = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(elife_provenance))]["UGRP_Lab_ID"].to_list()

### NOTE: these are samples for whom we extract FASTQ from CIDR CRAMs and realign. ###
realignment_provenance = (["CIDR_rd1", "CIDR_rd1,CIDR_rd2"])
samples_to_realign = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(realignment_provenance))]["UGRP_Lab_ID"].to_list()

### NOTE: these are samples for whom we align FASTQ extracted from ORA files, and treat the ### 
### resulting alignments as the final CRAMs. ###
alignment_provenance = (["UofU_rd1", "UofU_rd2,eLife", "UofU_rd2", "CIDR_rd1,UofU_rd2"])
samples_to_align = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(alignment_provenance))]["UGRP_Lab_ID"].to_list()

ad_hoc_provenance = (["CIDR_rd1,CIDR_rd2,UofU_rd1", "CIDR_rd1,CIDR_rd2,UofU_rd2"])
samples_ad_hoc = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(ad_hoc_provenance))]["UGRP_Lab_ID"].to_list()

print (samples_ad_hoc)

# find samples without sex info
no_sex = [s for s in samples_to_align if s not in SMP2SEX]

simple_samples = samples_to_keep + samples_to_realign + samples_to_align

samples = [s for s in simple_samples if s not in no_sex]
print (len(samples))

rule all:
    input:
        # expand("data/vcf/per-chrom/{SAMPLE}.GRCh38.{CHROM}.g.vcf.gz", SAMPLE = samples, CHROM=["chr1"]),
        expand("data/cram/{SAMPLE}.dupmarked.cram", SAMPLE = samples_ad_hoc)
        #expand("data/vcf/joint_genotyped/GRCh38.{CHROM}.joint_genotyped.bcf", CHROM=CHROMS)
        # expand("data/fastq/{SAMPLE}.{which}.fastq.gz", SAMPLE=samples_ad_hoc, which=["1", "2"])


rule joint_genotype:
    input:
        sif = "glnexus_v1.2.7.sif",
        gvcfs = expand("data/vcf/per-chrom/{SAMPLE}.{{ASSEMBLY}}.{{CHROM}}.g.vcf.gz", SAMPLE=samples),
    output: "data/vcf/joint_genotyped/{ASSEMBLY}.{CHROM}.joint_genotyped.bcf"
    threads: 16
    params:
        gl_nexus_prefix = lambda wildcards: f"gl_nexus_dbs/{wildcards.ASSEMBLY}_{wildcards.CHROM}"
    resources:
        mem_mb = 64_000
    script:
        "rules/bash_scripts/joint_genotype.sh"