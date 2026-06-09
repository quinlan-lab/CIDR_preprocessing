import pandas as pd
from collections import defaultdict
import glob
import json
import os

prov = pd.read_csv("PROVENANCE.tsv", sep="\t")

res = []
missing = []
for i, row in prov.iterrows():
    provenance = row["sequencing_provenance"].split(",")
    sample = row["UGRP_Lab_ID"]
    cram_path = None
    # if this sample is from eLife ONLY, we can leave the CRAM as-is
    if provenance == ["eLife"]:
        cram_path = f"/scratch/ucgd/lustre-core/UCGD_Datahub/Mosaic/920/UCGD/GRCh38/Data/PolishedCrams/{sample}.cram"
    # otherwise, we create a new CRAM
    else:
        cram_path = f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/cram/{sample}.cram"
    if cram_path is not None:
        if os.path.exists(cram_path):
            #if provenance != ["eLife"]:
            res.append({"ugrp_sample_id": sample, "cram_fh": cram_path})
        else:
            print (sample)
            missing.append({"ugrp_sample_id": sample, "provenance": ",".join(provenance)})

print(len(res))
print(pd.DataFrame(missing).groupby("provenance").size())

with open("json/cram_mapping.realigned.json", "w") as f:
    json.dump(res, f, indent=4)


### NOTE: here, we map samples to the CRAMs from which
# we'll need to extract FQ. these extracted FQ will be
# concatenated with FQ extracted from .ora files in the
# UofU_rd1/2 sequencing and aligned.

# read in sample information for the original CIDR bolus
CIDR_ORIG_INFO = (
    pd.read_excel(
        "Quinlan_Released_Data/Sample_Info/QuinlanNeklason_SIF.xlsx",
        sheet_name="Sheet0",
    ).dropna()
)[["Subject_ID", "Individual"]].rename(columns={"Subject_ID": "SUBJECT_ID"})

CIDR_ORIG_INFO["SUBJECT_ID"] = CIDR_ORIG_INFO["SUBJECT_ID"].astype(int).astype(str)

# read in mapping for original CIDR
CIDR_ORIG_MAP = pd.read_csv(
    "Quinlan_Released_Data/Sample_Info/SubjectSampleMappingFile_QuinlanNeklason.csv",
    dtype={"SUBJECT_ID": str},
)

CIDR_ORIG_INFO = CIDR_ORIG_INFO.merge(CIDR_ORIG_MAP, how="outer").rename(
    columns={
        "SUBJECT_ID": "UGRP_Lab_ID",
        "SAMPLE_ID": "prefix",
    }
)
CIDR_ORIG_INFO = CIDR_ORIG_INFO.dropna(subset=["prefix"])

# read in sample information for the topped-up CIDR bolus
CIDR_TOPUP_INFO = pd.read_csv(
    "dataset_to_PI_release2/samples_below_30x_with_generation_CIDRedit.csv",
    dtype={"Subject_ID": str},
)[
    [
        "Subject_ID",
        "wants more seq",
        "sample_id",
        "PICARD_average_alignment_coverage_over_genome",
    ]
].rename(
    columns={
        "Subject_ID": "UGRP_Lab_ID",
        "wants more seq": "topped_up",
        "sample_id": "prefix",
    }
)
# remove samples that didn't produce sequencing data
CIDR_TOPUP_INFO = CIDR_TOPUP_INFO[
    CIDR_TOPUP_INFO["PICARD_average_alignment_coverage_over_genome"]
    != "Library attempted but no sequence data generated"
]

res = []
for i, row in prov.iterrows():
    provenance = row["sequencing_provenance"].split(",")
    sample = row["UGRP_Lab_ID"]
    cram_path = None
    if provenance in (["CIDR_rd1"], ["CIDR_rd1", "UofU_rd2"]):
        prefix = CIDR_ORIG_INFO[CIDR_ORIG_INFO["UGRP_Lab_ID"] == sample][
            "prefix"
        ].to_list()
        assert len(prefix) == 1
        prefix = prefix.pop()
        cram_path = f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/Quinlan_Released_Data/CRAM/{prefix}.cram"

    elif provenance == ["CIDR_rd1", "CIDR_rd2"]:
        prefix = CIDR_TOPUP_INFO[CIDR_TOPUP_INFO["UGRP_Lab_ID"] == sample][
            "prefix"
        ].to_list()
        assert len(prefix) == 1
        prefix = prefix.pop()
        cram_path = f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/dataset_to_PI_release2/CRAM/{prefix}.cram"

    elif provenance in (
        ["CIDR_rd1", "CIDR_rd2", "UofU_rd1"],
        ["CIDR_rd1", "CIDR_rd2", "UofU_rd2"],
    ):
        prefix = CIDR_TOPUP_INFO[CIDR_TOPUP_INFO["UGRP_Lab_ID"] == sample][
            "prefix"
        ].to_list()
        assert len(prefix) == 1
        prefix = prefix.pop()
        cram_path = f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/dataset_to_PI_release2/CRAM/{prefix}.cram"

    elif provenance == ["UofU_rd2", "eLife"]:
        cram_path = f"/scratch/ucgd/lustre-core/UCGD_Datahub/Mosaic/920/UCGD/GRCh38/Data/PolishedCrams/{sample}.cram"

    #
    elif provenance == ["CIDR_rd1", "CIDR_rd1"]:
        cram_path = f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/Quinlan_Released_Data/CRAM/{sample}.cram"
    else:
        pass
#     if cram_path is not None and os.path.exists(cram_path):
#         res.append({"ugrp_sample_id": sample, "cram_fh": cram_path})
#     else:
#         missing.append({"ugrp_sample_id": sample, "provenance": ",".join(provenance)})
    if cram_path is not None:
        res.append({"ugrp_sample_id": sample, "cram_fh": cram_path})

# print (pd.DataFrame(missing).groupby("provenance").size())

with open("json/cram_mapping.to_extract.json", "w") as f:
    json.dump(res, f, indent=4)
