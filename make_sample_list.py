import pandas as pd

# get sample IDs in new CIDR CEPH cohort
sample_info = pd.read_csv(
    "Quinlan_Released_Data/Sample_Info/QuinlanNeklason_Master_Sample_Key.csv",
    dtype={
        "Subject_ID": "string",
    },
)

with open("sample_lists/CIDR_2025.samples.txt", "w") as outfh:
    for si, s in enumerate(sample_info["Subject_ID"].unique()):
        # if si > 5: break
        print (s.strip(), file=outfh)
    outfh.close()

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

with open("sample_lists/ELIFE_2019.samples.txt", "w") as outfh:
    for si, s in enumerate(orig_ped["SAMPLE_ID"].unique()):
        # if si not in (6, 9): continue
        print (s.strip(), file=outfh)
    outfh.close()
