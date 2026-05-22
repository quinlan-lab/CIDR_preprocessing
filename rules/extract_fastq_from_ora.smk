import pandas as pd
import glob

PREF = "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/"


### DEB ###
prefix2ora = {}
for fh in glob.glob(f"{PREF}UU_CIDR_CEPH/26866R/Fastq/*.ora"):
    prefix = fh.split("/")[-1].split("_")[0]
    prefix2ora[prefix] = fh

neklason_mapping = pd.read_excel(
    "data/2025 CEPH resequence submit core.xlsx",
    sheet_name="2025-08-29 CEPH REsequence list",
    dtype={"LINKID": str, "LABID": str},
)
# map each sample to their respective ora files
neklason_mapping["complete_file_path"] = neklason_mapping["26886Rx"].apply(lambda p: prefix2ora[f"26866X{p}"])
neklason_sample2ora = (
    neklason_mapping.groupby("LABID")
    .agg(ora_list=("complete_file_path", list))
    .to_dict()["ora_list"]
)

### SCOTT ###
watkins_mapping = pd.read_excel(
    f"{PREF}/UU_CIDR_CEPH/26501R_Id_list.xlsx",
    sheet_name="Sheet1",
)

manifest = pd.read_csv(f"{PREF}/UU_CIDR_CEPH/26501R/26501R_MANIFEST.csv")
manifest["Alt ID "] = manifest["sample_id"].apply(lambda s: "x" + s.split("X")[-1])

# merge
manifest = manifest.merge(watkins_mapping)
manifest["Sample Name"] = manifest["Sample Name"].astype(str)
manifest["complete_file_path"] = manifest["File"].apply(lambda fh: f"/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/UU_CIDR_CEPH/26501R/{fh}")

scott_sample2ora = manifest.groupby("Sample Name").agg(ora_list = ("complete_file_path", list)).to_dict()["ora_list"]

# combine scott and deb 
sample2ora = scott_sample2ora | neklason_sample2ora

print (sample2ora["150070"])

rule all:
    input:
        # expand("data/fastq/from_ora/{SAMPLE}.1.fastq.gz", SAMPLE = [k for k, v in sample2ora.items()])
        expand("data/fastq/from_ora/{SAMPLE}.1.fastq.gz", SAMPLE = ["110061", "110064"])

rule decompress:
    input:
        orad_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/orad.2.7.0.linux/orad",
        ora_list = lambda wildcards: [str(s) for s in sample2ora[wildcards.SAMPLE]],
    output:
        placeholder_fastq = "data/fastq/{SAMPLE}_0.fq.gz_1",
    params:
        n_input_files = lambda wildcards: len(sample2ora[wildcards.SAMPLE]) - 1,
        prefix = PREF
    threads: 16
    script:
        "bash_scripts/decompress_ora.sh"


rule concatenate:
    input:
        placeholder_fastq = "data/fastq/{SAMPLE}_0.fq.gz_1"
    output:
        fq1 = "data/fastq/from_ora/{SAMPLE}.1.fastq.gz",
        fq2 = "data/fastq/from_ora/{SAMPLE}.2.fastq.gz"
    params:
        n_input_files = lambda wildcards: len(sample2ora[wildcards.SAMPLE]),
        fastq1_list = lambda wildcards: [f"data/fastq/{wildcards.SAMPLE}_{n_fq}.fq.gz_1" for n_fq in range(len(sample2ora[wildcards.SAMPLE]))],
        fastq2_list = lambda wildcards: [f"data/fastq/{wildcards.SAMPLE}_{n_fq}.fq.gz_2" for n_fq in range(len(sample2ora[wildcards.SAMPLE]))]
    script:
        "bash_scripts/combine_fq.sh"