import json
import pandas as pd

CHROMS = [f"chr{c}" for c in range(1, 23)] + ["chrX"]

configfile: "config.yaml"

alignment_ref = config["alignment_ref"]
elife_cram_to_fastq_ref = config["elife_cram_to_fastq_ref"]
cidr_cram_to_fastq_ref = config["cidr_cram_to_fastq_ref"]
bam_to_cram_ref = config["bam_to_cram_ref"]

# we use a different cram2fastq reference for eLife and
# CIDR-aligned CRAMs.

neklason_mapping = pd.read_excel(
    "data/2025 CEPH resequence submit core.xlsx",
    sheet_name="2025-08-29 CEPH REsequence list",
    dtype={"LINKID": str, "LABID": str},
).rename(columns={"Unnamed: 0": "reseq_reason"})

SAMPLES_TO_MERGE = neklason_mapping[neklason_mapping["reseq_reason"].isin(["WASHU", "CIDRv2"])]["LABID"].unique()

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
with open("json/cram_mapping.realignedd.json") as f:
    dicts = json.load(f)
    for d in dicts:
        sample = d["ugrp_sample_id"]
        fh = d["cram_fh"]
        SMP2CRAM_NEW[sample] = fh

def get_new_cram_fh(wildcards):
    return SMP2CRAM_NEW[wildcards.SAMPLE]

def get_new_cram_idx_fh(wildcards):
    return SMP2CRAM_NEW[wildcards.SAMPLE] + ".crai"


wildcard_constraints:
    ASSEMBLY = "GRCh38",
    CHROM = r"chr[0-9]{1,2}|chrX|chrY",
    SAMPLE = r"[0-9]{4,7}"

elife_provenance = (["eLife"])
elife_samples = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(elife_provenance))]["UGRP_Lab_ID"].to_list()

SAMPLE2SOURCE = {}
SAMPLE2REF = {}

# NOTE: these are samples for whom we extract FASTQ from eLife CRAMs and
# merge with aligned data from ORA
elife_plus_provenance = (["UofU_rd2,eLife"])
_samples = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(elife_plus_provenance))]["UGRP_Lab_ID"].to_list()
for s in _samples:
    SAMPLE2SOURCE[s] = "merged"
    # which CRAM reference are we using to extract FASTQ?
    SAMPLE2REF[s] = elife_cram_to_fastq_ref

# NOTE: these are samples for whom we extract FASTQ from CIDR CRAMs 
# and merge with aligned data from ORA. CIDR rd2 superceedes rd1 if available,
# because rd2 is merged with rd1.
cidr_plus_provenance = (["CIDR_rd1,UofU_rd2", "CIDR_rd1,CIDR_rd2,UofU_rd1", "CIDR_rd1,CIDR_rd2,UofU_rd2"])
_samples = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(cidr_plus_provenance))]["UGRP_Lab_ID"].to_list()
for s in _samples:
    SAMPLE2SOURCE[s] = "merged"
    # which CRAM reference are we using to extract FASTQ?
    SAMPLE2REF[s] = cidr_cram_to_fastq_ref

# NOTE: these are samples for whom we extract FASTQ from CIDR CRAMs
# and simply align to a common reference.
realignment_provenance = (["CIDR_rd1", "CIDR_rd1,CIDR_rd2", "CIDR_rd1,CIDR_rd1"])
_samples = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(realignment_provenance))]["UGRP_Lab_ID"].to_list()
for s in _samples:
    SAMPLE2SOURCE[s] = "from_cram"
    # which CRAM reference are we using to extract FASTQ?
    SAMPLE2REF[s] = cidr_cram_to_fastq_ref

# NOTE: these are samples for whom we align FASTQ extracted from ORA files.
alignment_provenance = (["UofU_rd1", "UofU_rd2"])
_samples = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(alignment_provenance))]["UGRP_Lab_ID"].to_list()
for s in _samples:
    SAMPLE2SOURCE[s] = "from_ora"

def get_cram2fastq_ref(wildcards):
    return SAMPLE2REF[wildcards.SAMPLE]

include: "rules/merge_uurd2_with_ceph.smk"
include: "rules/cram2fastq.smk"
include: "rules/fastq2bam.smk"
include: "rules/call_variants.smk"

# samples for testing relatedness
som_testing = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(["eLife"]))]["UGRP_Lab_ID"].to_list()[:25]
som_testing.extend(["10092", "390070"])

cidr_samples = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(["CIDR_rd1", "CIDR_rd1,CIDR_rd2", "CIDR_rd1,CIDR_rd1"]))]["UGRP_Lab_ID"].to_list()
cidr_uofu_samples = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(["CIDR_rd1,UofU_rd2", "CIDR_rd1,CIDR_rd2,UofU_rd1", "CIDR_rd1,CIDR_rd2,UofU_rd2"]))]["UGRP_Lab_ID"].to_list()
elife_uofu_samples = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(["UofU_rd2,eLife"]))]["UGRP_Lab_ID"].to_list()
uofu_samples = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(["UofU_rd1", "UofU_rd2"]))]["UGRP_Lab_ID"].to_list()
# samples = PROVENANCE[(PROVENANCE["sequencing_provenance"].isin(["CIDR_rd1"]))]["UGRP_Lab_ID"].to_list()
# all_samples = [s for s in all_samples if s  not in ("130053", "210071", "300010")]
samples = uofu_samples + cidr_uofu_samples + elife_uofu_samples + cidr_samples
# samples = [s for s in samples if s in SMP2SEX]
samples = [s for s in samples if s != "NA12878"]
rule all:
    input:
        # expand("data/vcf/per-chrom/{SAMPLE}.GRCh38.{CHROM}.g.vcf.gz", SAMPLE = ["150070", "400050", "80010", "300100"], CHROM=CHROMS),
        expand("/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/cram/{SAMPLE}.cram.crai", SAMPLE = samples),
        # expand("data/fastq/merged/{SAMPLE}.1.fastq.gz", SAMPLE = ["8321"])
        # expand("data/vcf/joint_genotyped/GRCh38.{CHROM}.joint_genotyped.bcf", CHROM=["chr1"])
        # expand("data/fastq/{SAMPLE}.{which}.fastq.gz", SAMPLE=samples_ad_hoc, which=["1", "2"])
        # "data/vcf/merged/GRCh38.merged.bcf"
        # "somalier.html"

rule run_somalier:
    input:
        cram = get_new_cram_fh,
        ref = bam_to_cram_ref
    output: "somalier_data/{SAMPLE}.somalier"
    shell:
        """
        somalier extract -d somalier_data/ \
                         --sites sites.hg38.vcf.gz \
                         -f {input.ref} \
                         {input.cram}
        """

rule relate_somalier:
    input:
        fhs = expand("somalier_data/{SAMPLE}.somalier", SAMPLE=som_testing),
        ped = "data/ceph.ped"
    output:
        "somalier.html"
    shell:
        """
        somalier relate --ped {input.ped} {input.fhs}
        """

rule joint_genotype:
    input:
        sif = "glnexus_v1.2.7.sif",
        gvcfs = expand("data/vcf/per-chrom/{SAMPLE}.{{ASSEMBLY}}.{{CHROM}}.g.vcf.gz", SAMPLE=samples),
    output: "data/vcf/joint_genotyped/{ASSEMBLY}.{CHROM}.joint_genotyped.bcf"
    threads: 16
    params:
        gl_nexus_prefix = lambda wildcards: f"gl_nexus_dbs/{wildcards.ASSEMBLY}_{wildcards.CHROM}"
    resources:
        mem_mb = 64_000,
        slurm_account = "ucgd-rw",
        slurm_partition = "ucgd-rw"
    script:
        "rules/bash_scripts/joint_genotype.sh"


rule combine_vcfs:
    input:
        vcfs = expand("data/vcf/joint_genotyped/{{ASSEMBLY}}.{CHROM}.joint_genotyped.bcf", CHROM=CHROMS)
    output:
        vcf = "data/vcf/merged/{ASSEMBLY}.merged.bcf"
    threads: 16
    resources:
        runtime = 1440
    shell:
        """
        module load bcftools
        
        bcftools concat --threads {threads} \
                        -Ob \
                        -o {output.vcf} \
                        {input.vcfs}
        """