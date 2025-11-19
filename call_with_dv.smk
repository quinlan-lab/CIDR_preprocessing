import json
import pandas as pd


CHROMS = list(map(str, range(1, 23)))
CHROMS = [f"chr{c}" for c in CHROMS]
CHROMS.extend(["chrX", "chrY"])

SMP2CRAM = {}
# map samples to CRAMs
with open("json/cram_mapping.json") as f:
    dicts = json.load(f)
    for d in dicts:
        sample = d["sample"]
        fh = d["cram_fh"]
        SMP2CRAM[sample] = fh
        
# get sample IDs in new CIDR CEPH cohort
sample_info = pd.read_csv(
    "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/ped/joint.cidr.ceph.ped",
    sep="\t",
    dtype={"SAMPLE_ID": str}
)
CIDR_SAMPLE_MAP = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/Quinlan_Released_Data/Sample_Info/SubjectSampleMappingFile_QuinlanNeklason.csv"
cidr_samples = pd.read_csv(CIDR_SAMPLE_MAP, dtype={"SUBJECT_ID": str})

# get 10 random trios
kids = sample_info[(sample_info["FATHER_ID"] != 0) & (sample_info["MOTHER_ID"] != 0)]
cidr_kids = kids[kids["SAMPLE_ID"].isin(cidr_samples["SUBJECT_ID"])]

cidr_trios = []
for i,row in cidr_kids.iterrows():
    for _id in ("SAMPLE_ID", "FATHER_ID", "MOTHER_ID"):
        cidr_trios.append(row[_id])

smp2sex = dict(zip(sample_info["SAMPLE_ID"], sample_info["SEX"]))

orig_ped = "ped/16-08-06_WashU-Yandell-CEPH.ped"

# read in the Ped file from the original CEPH 2019 study
orig_ped = pd.read_csv(
    orig_ped,
    sep="\t",
    names=[
        "FAMILY_ID",
        "SAMPLE_ID",
        "FATHER_ID",
        "MOTHER_ID",
        "SEX",
        "PHENOTYPE",
        "MISC",
    ],
    dtype={"FAMILY_ID": "string", "SAMPLE_ID": "string"},
)
# elife_smp2sex = dict(zip(orig_ped["SAMPLE_ID"], orig_ped["SEX"]))

# smp2sex.update(elife_smp2sex)

elife_kids = orig_ped[(orig_ped["FATHER_ID"] != 0) & (orig_ped["MOTHER_ID"] != 0)]

elife_trios = []
for i,row in elife_kids.iterrows():
    for _id in ("SAMPLE_ID", "FATHER_ID", "MOTHER_ID"):
        elife_trios.append(row[_id])

all_samples = list(map(str, elife_trios + cidr_trios))
all_samples = [s for s in all_samples if s in SMP2CRAM and s in smp2sex]
for s in all_samples:
    if str(s) not in SMP2CRAM:
        print (s)

print (len(set(all_samples)))

wildcard_constraints:
    ASSEMBLY = "GRCh38",
    CHROM = r"chr[0-9]{1,2}|chrX|chrY"

rule all:
    input:
        # expand("data/vcf/per-chrom/{SAMPLE}.GRCh38.{CHROM}.g.vcf.gz", SAMPLE=list(set(all_samples)), CHROM=CHROMS)
        "data/vcf/merged/GRCh38.joint_genotyped.vcf.gz.tbi",
        # "data/vcf/per-chrom/GRCh38.chr21.joint_genotyped.vcf.gz"

rule call_snvs:
    input:
        ref = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        ref_idx = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz.fai",
        cram = lambda wildcards: SMP2CRAM[wildcards.SAMPLE],
        sif = "deepvariant_1.9.0.sif"
    output:
        vcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf.gz",
        gvcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.g.vcf.gz",
    threads: 8
    params:
        haploid_contigs_arg = lambda wildcards: "--haploid_contigs chrX,chrY" if smp2sex[wildcards.SAMPLE] == 1 and wildcards.CHROM in ("chrX", "chrY") else ""
    resources:
        mem_mb = 32_000,
    script:
        "bash_scripts/run_deepvariant.sh"


rule combine_trio_chrom_vcfs:
    input:
        sif = "glnexus_v1.2.7.sif",
        gvcfs = expand("data/vcf/per-chrom/{SAMPLE}.{{ASSEMBLY}}.{{CHROM}}.g.vcf.gz", SAMPLE=list(set(all_samples))),
    output: "data/vcf/per-chrom/{ASSEMBLY}.{CHROM}.joint_genotyped.vcf.gz"
    threads: 8
    params:
        gl_nexus_prefix = lambda wildcards: f"gl_nexus_dbs/{wildcards.ASSEMBLY}_{wildcards.CHROM}"
    resources:
        mem_mb = 64_000
    script:
        "bash_scripts/joint_genotype.sh"


rule merge_trio_vcfs:
    input:
        vcfs = expand("data/vcf/per-chrom/{{ASSEMBLY}}.{CHROM}.joint_genotyped.vcf.gz", CHROM=CHROMS)
    output: "data/vcf/merged/{ASSEMBLY}.joint_genotyped.vcf.gz"
    threads: 16
    resources:
        runtime = 1440
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} | bgzip > {output}
        """


rule index_merged_vcf:
    input: vcf = "data/vcf/merged/{ASSEMBLY}.joint_genotyped.vcf.gz"
    output: "data/vcf/merged/{ASSEMBLY}.joint_genotyped.vcf.gz.tbi"
    script:
        "bash_scripts/index_vcf.sh"