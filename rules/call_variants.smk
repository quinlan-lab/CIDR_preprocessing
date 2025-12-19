wildcard_constraints:
    ASSEMBLY = "GRCh38",
    CHROM = r"chr[0-9]{1,2}|chrX|chrY"


ELIFE_REF_FH = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"

# get sample IDs in new CIDR CEPH cohort
sample_info = pd.read_csv(
    "/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/ped/joint.cidr.ceph.ped",
    sep="\t",
    dtype={"SAMPLE_ID": str}
)
SMP2SEX = dict(zip(sample_info["SAMPLE_ID"], sample_info["SEX"]))


# map sample names to CRAM files
SMP2CRAM = {}
with open("json/cram_mapping.json") as f:
    dicts = json.load(f)
    for d in dicts:
        sample = d["sample"]
        fh = d["cram_fh"]
        SMP2CRAM[sample] = fh


rule call_snvs:
    input:
        ref = ELIFE_REF_FH,
        cram = lambda wildcards: SMP2CRAM[wildcards.SAMPLE],
        sif = "deepvariant_1.9.0.sif"
    output:
        vcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.vcf.gz",
        gvcf = "data/vcf/per-chrom/{SAMPLE}.{ASSEMBLY}.{CHROM}.g.vcf.gz",
    threads: 8
    params:
        haploid_contigs_arg = lambda wildcards: "--haploid_contigs chrX,chrY" if SMP2SEX[wildcards.SAMPLE] == 1 and wildcards.CHROM in ("chrX", "chrY") else ""
    resources:
        mem_mb = 32_000,
    script:
        "bash_scripts/run_deepvariant.sh"