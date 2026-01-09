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


rule all:
    input:
        expand("data/vcf/per-sample/{SAMPLE}.GRCh38.g.vcf.gz", SAMPLE = elife),
        expand("data/cram/{SAMPLE}.dupmarked.cram", SAMPLE = cidr)