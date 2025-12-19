include: "rules/cram2fastq.smk"
include: "rules/fastq2bam.smk"

import json
import pandas as pd


CIDR_SAMPLE_MAP = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/Quinlan_Released_Data/Sample_Info/SubjectSampleMappingFile_QuinlanNeklason.csv"
cidr_samples = pd.read_csv(CIDR_SAMPLE_MAP, dtype={"SUBJECT_ID": str})


rule all:
    input:
        "data/cram/230091.dupmarked.cram.crai"