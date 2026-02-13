import pandas as pd
from collections import defaultdict
import glob
import json
import os

PREF = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH"
MASTER_PED = "/scratch/ucgd/lustre-labs/quinlan/u0890814/CIDR_4Gen/ped_files/master_ped_all_info.ped"

# get sample info for full cohort from Julia's "ped"
SAMPLE_INFO = pd.read_csv(
    MASTER_PED,
    sep="\t",
    dtype={"UGRP_Lab_ID": str, "Gender": str}
)
SAMPLE_INFO = SAMPLE_INFO[SAMPLE_INFO["sequencing_available"] == "yes"]
print (SAMPLE_INFO.shape)

# print (SAMPLE_INFO[SAMPLE_INFO["UGRP_Lab_ID"] == "200095"][["Sequencing"]])

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

# read in sample information for the topped-up CIDR bolus
CIDR_TOPUP_INFO = pd.read_csv(
    "dataset_to_PI_release2/samples_below_30x_with_generation_CIDRedit.csv",
    dtype={"Subject_ID": str},
)[
    [
        "Subject_ID",
        "wants more seq",
        "sample_id",
    ]
].rename(
    columns={
        "Subject_ID": "UGRP_Lab_ID",
        "wants more seq": "topped_up",
        "sample_id": "prefix",
    }
)

# merge the original CIDR bolus with the topped-up CIDR bolus so we know who was topped up
CIDR_MERGED = CIDR_ORIG_INFO.merge(CIDR_TOPUP_INFO, how="left").fillna("no")
CIDR_MERGED.replace(to_replace={"wants more seq": "yes", "no": "NA"}, inplace=True)

# NOTE: six of the samples in the CIDR topup bolus had no sequence data generated
# for sample in CIDR_MERGED.query("topped_up == 'yes'")["topup_prefix"].unique():
#     if not os.path.exists(f"dataset_to_PI_release2/CRAM/{sample}.cram"):
#         print (sample, CIDR_MERGED[CIDR_MERGED["topup_prefix"] == sample])

# read in Deb's bolus of sequencing metadata
DEB_INFO = pd.read_excel(
    "data/2025 CEPH resequence submit core.xlsx",
    sheet_name="2025-08-29 CEPH REsequence list",
    dtype={"LABID": str},
)[["LABID"]].rename(columns={"LABID": "UGRP_Lab_ID"})

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

ELIFE_INFO = SAMPLE_INFO[SAMPLE_INFO["Sequencing"] == "WashU-Illumina_short-read"][["UGRP_Lab_ID"]]

# for every sample in the master ped, figure out whether we sequenced it in
# a) the original eLife
# b) the first CIDR bolus
# c) the second CIDR bolus
# d) scott's bolus
# e) deb's bolus
# then, figure out which CRAM file already exists for this sample (if any) and
# whether we need to realign it.

res = []
for i, row in SAMPLE_INFO.iterrows():
    sample_id = row["UGRP_Lab_ID"]

    # figure out the membership of this sample in each bolus
    in_cidr = sample_id in CIDR_MERGED["UGRP_Lab_ID"].to_list()
    in_elife = sample_id in ELIFE_INFO["UGRP_Lab_ID"].to_list()
    in_deb = sample_id in DEB_INFO["UGRP_Lab_ID"].to_list()
    in_scott = sample_id in WATKINS_INFO["UGRP_Lab_ID"].to_list()

    # how many sequencing runs was this sample in?
    overlap = in_cidr + in_elife + in_deb + in_scott
    assert 1 <= overlap <= 2

    # if this sample is only seen in one dataset, it's relatively easy to determine what
    # the "final" CRAM path should be. we always realign w/r/t the original eLife samples,
    # so if this sample is only present in the original eLife dataset, we'll just use their
    # original CRAM
    if overlap == 1:
        if in_elife:
            final_cram_path = f"/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/CEPH/cram/{sample_id}.cram"
            realign = False
            fastq_path = ""
            orig_cram_path = ""
        elif in_cidr:
            final_cram_path = f"{PREF}/data/cram/{sample_id}.dupmarked.cram"
            realign = True
            fastq_path = f"{PREF}/data/fastq/{sample_id}.1.fastq.gz,{PREF}/data/fastq/{sample_id}.2.fastq.gz"

            cidr_entry = CIDR_MERGED[CIDR_MERGED["UGRP_Lab_ID"] == sample_id]            
            prefixes = cidr_entry["prefix"].to_list()
            topped_up = cidr_entry["topped_up"].to_list()
            if sample_id in ("300100", "80010"):
                assert len(prefixes) == 2
                prefixes = [sample_id]
            assert len(set(topped_up)) == 1
            if topped_up.pop() == "yes":
                orig_cram_path = ",".join([f"{PREF}/dataset_to_PI_release2/CRAM/{prefix}.cram" for prefix in prefixes])
            else:
                orig_cram_path = ",".join([f"{PREF}/Quinlan_Released_Data/CRAM/{prefix}.cram" for prefix in prefixes])
        else:
            final_cram_path = f"{PREF}/data/cram/{sample_id}.dupmarked.cram"
            realign = True
            fastq_path = f"{PREF}/data/fastq/{sample_id}.1.fastq.gz,{PREF}/data/fastq/{sample_id}.2.fastq.gz"
            orig_cram_path = ""

    else:
        orig_cram_path = ""
        if in_elife:
            orig_cram_path = f"/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/CEPH/cram/{sample_id}.cram"
        elif in_cidr:
            cidr_entry = CIDR_MERGED[CIDR_MERGED["UGRP_Lab_ID"] == sample_id]
            prefixes = cidr_entry["prefix"].to_list()
            topped_up = cidr_entry["topped_up"].to_list()
            assert len(set(topped_up)) == 1
            if topped_up.pop() == "yes":
                orig_cram_path = ",".join([f"{PREF}/dataset_to_PI_release2/CRAM/{prefix}.cram" for prefix in prefixes])
            else:
                orig_cram_path = ""
        final_cram_path = f"{PREF}/data/cram/{sample_id}.dupmarked.cram"
        realign = True
        fastq_path = f"{PREF}/data/fastq/{sample_id}.1.fastq.gz,{PREF}/data/fastq/{sample_id}.2.fastq.gz"

    res.append({
        "UGRP_Lab_ID": sample_id,
        "REALIGN": realign,
        "in_elife": in_elife,
        "in_cidr": in_cidr,
        "in_gnomex": in_deb,
        "in_watkins": in_scott,
        "final_cram_path": final_cram_path,
        "orig_cram_path": orig_cram_path,
        "fastq_paths": fastq_path,
    })

res_df = pd.DataFrame(res).replace({"": "NA"})

to_realign = res_df[res_df["REALIGN"] == True]
# make sure that samples in two datasets have an original cram path
# print (to_realign.query("overlap == 2 and orig_cram_path == 'NA'"))
res_df.to_csv("o.tsv", sep="\t", index=False)
# # what about the "original" CRAM paths from which we'll extract FASTQ?
# if overlap == 1 and in_elife:

# elif in_cidr:
#     cidr_entry = CIDR_MERGED[CIDR_MERGED["UGRP_Lab_ID"] == sample_id]
#     prefixes = cidr_entry["prefix"].to_list()
#     topped_up = cidr_entry["topped_up"].to_list()
#     assert len(set(topped_up)) == 1
#     if topped_up.pop() == "yes":
#         cram_path = ",".join([f"{PREF}/dataset_to_PI_release2/CRAM/{prefix}.cram" for prefix in prefixes])
#     else:
#         cram_path = ",".join([f"{PREF}/Quinlan_Released_Data/CRAM/{prefix}.cram" for prefix in prefixes])
#     realign = True
# # none of Deb's samples should be unique to the sequencing bolus.
# # they should all be resequencing/topup of a sample with an existing CRAM
# else:
#     pass
# if there's overlap, we'll need to align FASTQs "from scratch," so there's
# no existing CRAM path
# else:
#     cram_path = ""
#     # if this sample is an original eLife sample, but was sequenced in the
#     # new deb bolus, we need to realign
#     if in_elife and in_deb:
#         cram_path = ""
#     # if this sample was in CIDR and sequenced by Deb, deb's sequencing takes precedent
#     elif in_cidr and in_deb:
#         cram_path = ""
#     # if this sample was in CIDR and sequence by scott, scott takes precedent ??? NOTE
#     elif in_cidr and in_scott:
#         cram_path = ""
#     else:
#         print ("!! AH !!")
# return cram_path, realign


# # if we're in the CIDR boli, figure out which
# # of the two CIDR boli should supercede
# if sample_id in CIDR_MERGED["UGRP_Lab_ID"]:
#     cidr_entry = CIDR_MERGED[CIDR_MERGED["UGRP_Lab_ID"] == sample_id]
#     assert cidr_entry.shape[0] == 1
#     prefix = cidr_entry["prefix"].to_list().pop()
#     topped_up = cidr_entry["topped_up"].to_list().pop()
#     if topped_up == "yes":
#         cram_path = f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/dataset_to_PI_release2/CRAM/{prefix}.cram"
#         realign = True
#     else:
#         cram_path = f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/Quinlan_Released_Data/CRAM/{prefix}.cram"
#         realign = True
# # if we're in scott's bolus
# # otherwise, if this is an existing eLife sample
# # if we're in the deb bolus, we need to check if we're going to supercede
# # any eLife samples
# elif sample_id in DEB_INFO["UGRP_Lab_ID"]:
#     if sample_id in ELIFE_INFO["UGRP_Lab_ID"]:

# SAMPLE_INFO["original_cram_path"] = SAMPLE_INFO["UGRP_Lab_ID"].apply(lambda s: get_original_cram_paths(s)[0])
# print ("#")
# SAMPLE_INFO["realign"] = SAMPLE_INFO["UGRP_Lab_ID"].apply(lambda s: get_original_cram_paths(s)[1])
# print (SAMPLE_INFO.groupby(["original_cram_path", "realign"]).size())
# PREF = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/"
# # path to 2019 CEPH data
# ORIG_PREF = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/CEPH"

# out = []

# # get the metadata for deb's bolus of sequencing data
# neklason_metadata = pd.read_excel(
#     "data/2025 CEPH resequence submit core.xlsx",
#     sheet_name="2025-08-29 CEPH REsequence list",
#     dtype={"LABID": str},
# )
# # some of the original eLife samples were resequenced, and we'll want to
# # use their updated CRAMs for variant calling
# for sample in SAMPLE_INFO["UGRP_Lab_ID"].unique():

#     # we have existing CRAM files for the 603 original eLife samples.
#     # everyone else was resequenced.
#     cram_for_gvcf = f"{ORIG_PREF}/cram/{sample}.cram"
#     if sample in neklason_metadata["LABID"].unique():
#         cram_for_gvcf = f"data/cram/{sample}.dupmarked.cram"
#     out.append({"UGRP_LABID": sample, "cram_for_gvcf": cram_for_gvcf})

# assert len(out) == 603

# print (SAMPLE_INFO.groupby("Sequencing").size())
# # output file with mappings from sample names to GVCFs
# sample_to_cram = []

# watkins_mapping = pd.read_excel(f"{PREF}/UU_CIDR_CEPH/26501R_Id_list.xlsx", sheet_name="Sheet1")
# manifest = pd.read_csv(f"{PREF}/UU_CIDR_CEPH/26501R/26501R_MANIFEST.csv")
# manifest["Alt ID "] = manifest["sample_id"].apply(lambda s: "x" + s.split("X")[-1])

# # merge
# manifest = manifest.merge(watkins_mapping)
# manifest["Sample Name"] = manifest["Sample Name"].astype(str)
# print (manifest[manifest["Sample Name"] == "210071"])
# print (manifest["Sample Name"].nunique())
# for s in manifest["Sample Name"].unique():
#     cram_fh = f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/cram/{s}.dupmarked.cram"
#     if not os.path.exists(cram_fh):
#         print (s)
#         continue
#     sample_to_cram.append({"ugrp_sample_id": s, "cram_fh": f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/cram/{s}.dupmarked.cram"})

# MASTER_PED = "/scratch/ucgd/lustre-labs/quinlan/u0890814/CIDR_4Gen/ped_files/master_ped_all_info.ped"

# # get sample IDs in new CIDR CEPH cohort
# sample_info = pd.read_csv(
#     MASTER_PED,
#     sep="\t",
#     dtype={"UGRP_Lab_ID": str, "Gender": int}
# )

# orig = sample_info[sample_info["Sequencing"] == "WashU-Illumina_short-read"]

# orig_sample_ids = orig["UGRP_Lab_ID"].unique()

# # create mapping of original CEPH sample IDs to GVCF file names
# smp2fh = defaultdict()
# for fh in glob.glob(f"{ORIG_PREF}/cram/*.cram"):
#     smp = fh.split("/")[-1].split(".")[0]
#     for s in orig_sample_ids:
#         if s == smp:
#             sample_to_cram.append({"ugrp_sample_id": s, "cram_fh": fh})


# cidr = sample_info[sample_info["Sequencing"].isin(["CIDR-Illumina_short-read", "CIDR-Illumina_short-read_top-up"])]
# # create a dictionary
# for i, row in cidr.iterrows():
#     cram_fh = f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/data/cram/{row['UGRP_Lab_ID']}.dupmarked.cram"
#     if not os.path.exists(cram_fh):
#         print (row["UGRP_Lab_ID"])
#     sample_to_cram.append({"ugrp_sample_id": row["UGRP_Lab_ID"], "cram_fh": cram_fh})

# print (len(sample_to_cram))

# with open("json/cram_mapping.realigned.json", "w") as f:
#     json.dump(sample_to_cram, f, indent=4)
