import pandas as pd

PREF = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH"
MASTER_PED = "/scratch/ucgd/lustre-labs/quinlan/u0890814/CIDR_4Gen/ped_files/master_ped_all_info.ped"

# get sample info for full cohort from Julia's "ped"
SAMPLE_INFO = pd.read_csv(
    MASTER_PED,
    sep="\t",
    dtype={"UGRP_Lab_ID": str, "Gender": str}
)

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
    dtype={"SUBJECT_ID": str}
)

CIDR_ORIG_INFO = CIDR_ORIG_INFO.merge(CIDR_ORIG_MAP, how="outer").rename(
    columns={
        "SUBJECT_ID": "UGRP_Lab_ID",
        "SAMPLE_ID": "prefix",
    }
)
CIDR_ORIG_INFO = CIDR_ORIG_INFO.dropna(subset=["prefix"])

CIDR_ORIG_INFO["provenance"] = "CIDR_rd1"

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
CIDR_TOPUP_INFO = CIDR_TOPUP_INFO[CIDR_TOPUP_INFO["PICARD_average_alignment_coverage_over_genome"] != "Library attempted but no sequence data generated"]

CIDR_TOPUP_INFO["provenance"] = "CIDR_rd2"


# read in Deb's bolus of sequencing metadata
DEB_INFO = pd.read_excel(
    "data/2025 CEPH resequence submit core.xlsx",
    sheet_name="2025-08-29 CEPH REsequence list",
    dtype={"LABID": str},
)[["LABID"]].rename(columns={"LABID": "UGRP_Lab_ID", "Unnamed: 0": "reseq_reason"})

DEB_INFO["provenance"] = "UofU_rd2"

# read in scott's bolus of sequencing metadata
WATKINS_INFO = pd.read_excel(
    "UU_CIDR_CEPH/26501R_Id_list.xlsx",
    sheet_name="Sheet1",
    dtype={"Sample Name": str},
).rename(
    columns={
        "Sample Name": "UGRP_Lab_ID",
        "Alt ID ": "topup_prefix",
    }
)[
    ["UGRP_Lab_ID", "topup_prefix"]
]
WATKINS_INFO["provenance"] = "UofU_rd1"


ELIFE_INFO = SAMPLE_INFO[SAMPLE_INFO["Sequencing"] == "WashU-Illumina_short-read"][["UGRP_Lab_ID"]]
ELIFE_INFO["provenance"] = "eLife"

merged = pd.concat([CIDR_ORIG_INFO, CIDR_TOPUP_INFO, DEB_INFO, WATKINS_INFO, ELIFE_INFO])
prov = merged.groupby("UGRP_Lab_ID").agg(sequencing_provenance = ("provenance", lambda p: ",".join(p)), n_prov = ("provenance", lambda p: len(p))).reset_index()
prov.groupby("sequencing_provenance").size().reset_index().rename(columns={0: "count"}).sort_values("count", ascending=False).to_csv("a.tsv", sep="\t", index=False)

# these samples should be treated as being only UofUrd2, since we don't want to merge
CONTAMINATED_WASHU = """
8526
1416
8092
1360
8510
8530
8504
8612
8512
8524""".split()

prov["sequencing_provenance"] = prov.apply(
    lambda row: (
        "UofU_rd2"
        if row["UGRP_Lab_ID"] in CONTAMINATED_WASHU
        else row["sequencing_provenance"]
    ),
    axis=1,
)
prov.to_csv("PROVENANCE.tsv", sep="\t", index=False)
