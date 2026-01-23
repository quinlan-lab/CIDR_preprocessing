include: "rules/call_variants.smk"
include: "rules/cram2fastq.smk"
include: "rules/fastq2bam.smk"

import json
import pandas as pd

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

# cidr_dads = cidr["Father"].to_list()
# cidr_moms = cidr["Mother"].to_list()
# cidr_kids = cidr["UGRP_Lab_ID"].to_list()

# samples = cidr_dads + cidr_moms + cidr_kids
# print (len(set(cidr_kids)))
# cidr_kids = ["200092"]

CHROMS = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]

rule all:
    input:
        expand("data/vcf/joint_genotyped/{ASSEMBLY}.{CHROM}.joint_genotyped.vcf.gz", ASSEMBLY=["GRCh38"], CHROM=CHROMS)
        # expand("data/vcf/per-sample/{SAMPLE}.GRCh38.g.vcf.gz", SAMPLE = list(set(cidr_kids))[:10]),
        # expand("data/cram/{SAMPLE}.dupmarked.cram.crai", SAMPLE = cidr)


rule joint_genotype:
    input:
        sif = "glnexus_v1.2.7.sif",
        gvcfs = expand("data/vcf/per-chrom/{SAMPLE}.{{ASSEMBLY}}.{{CHROM}}.g.vcf.gz", SAMPLE=list(set(cidr + elife))),
    output: "data/vcf/joint_genotyped/{ASSEMBLY}.{CHROM}.joint_genotyped.vcf.gz"
    threads: 16
    params:
        gl_nexus_prefix = lambda wildcards: f"gl_nexus_dbs/{wildcards.ASSEMBLY}_{wildcards.CHROM}"
    resources:
        mem_mb = 64_000
    script:
        "rules/bash_scripts/joint_genotype.sh"


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