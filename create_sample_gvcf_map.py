import pandas as pd
from collections import defaultdict
import glob
import json


# output file with mappings from sample names to GVCFs
sample_to_gvcf = []

# path to CIDR GVCFs
CIDR_GVCF_PREF = "Quinlan_Released_Data/GVCF"
# path to 2019 CEPH data
ORIG_PREF = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/CEPH"

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
orig_sample_ids = orig_ped["SAMPLE_ID"].unique()

# create mapping of original CEPH sample IDs to GVCF file names
smp2fh = defaultdict()
for fh in glob.glob(f"{ORIG_PREF}/vcf/snv_indels/*.gz"):
    for s in orig_sample_ids:
        if s in fh:
            sample_to_gvcf.append({"sample": s, "gvcf_fh": fh})


# get sample IDs in new CIDR CEPH cohort
sample_info = pd.read_csv(
    "Quinlan_Released_Data/Sample_Info/QuinlanNeklason_Master_Sample_Key.csv",
    dtype={
        "Subject_ID": "string",
    },
)

# create a dictionary 
for fh in glob.glob(f"{CIDR_GVCF_PREF}/*.gz"):
    for i, row in sample_info.iterrows():
        if row["SM_Tag"] in fh:
            sample_to_gvcf.append({"sample": row["Subject_ID"], "gvcf_fh": fh})


with open("json/gvcf_mapping.json", "w") as f:
    json.dump(sample_to_gvcf, f, indent=4)
