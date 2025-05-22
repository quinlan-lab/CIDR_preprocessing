import pandas as pd
import glob
from collections import defaultdict
import numpy as np

CIDR_GVCF_PREF = "Quinlan_Released_Data/GVCF"
ORIG_PREF = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/CEPH"
FAMILY_ID = "1463"
DB = "/scratch/ucgd/lustre-labs/quinlan/u1006375/my_database/"
REF = "/scratch/ucgd/lustre/common/data/Reference/homo_sapiens/GRCh38/primary_assembly_decoy_phix_masked.fa"

CHROMS = list(map(str, range(1, 23)))
CHROMS = [f"chr{c}" for c in CHROMS]


orig_ped = f"{ORIG_PREF}/ped/2019-6-27-manual-checked-ceph.no_missing_samples.sorted.peddy_filtered.ped"

# get sample IDs in original CEPH cohort
orig_ped = pd.read_csv(orig_ped, sep="\t", names=["FAMILY_ID", "SAMPLE_ID", "FATHER_ID", "MOTHER_ID", "SEX", "PHENOTYPE", "AGE", "DENOMINATOR"], dtype={"FAMILY_ID": "string", "SAMPLE_ID": "string"})
orig_ped = orig_ped[orig_ped["FAMILY_ID"] == FAMILY_ID]
orig_sample_ids = orig_ped["SAMPLE_ID"].unique()

# create mapping of sample IDs to file names
smp2fh = defaultdict()
for fh in glob.glob(f"{ORIG_PREF}/vcf/snv_indels/*.gz"):
    for s in orig_sample_ids:
        if s in fh:
          smp2fh[s] = fh

# get sample IDs in new CIDR CEPH cohort
sample_info = pd.read_csv("Quinlan_Released_Data/Sample_Info/SubjectSampleMappingFile_QuinlanNeklason.csv", dtype={"SUBJECT_ID": "string", "SAMPLE_ID": "string"},)
cidr_ped = pd.read_excel("Quinlan_Released_Data/Sample_Info/QuinlanNeklason_SIF.xlsx", sheet_name="Sheet0", dtype={"Family": "string", "Individual": "string"})
cidr_ped = cidr_ped[cidr_ped["Family"] == FAMILY_ID].dropna()

cidr_sample_ids = cidr_ped["Individual"].unique()

sample_info = sample_info[sample_info["SUBJECT_ID"].isin(cidr_sample_ids)]
subject2sample = dict(zip(sample_info["SUBJECT_ID"], sample_info["SAMPLE_ID"]))
sample2subject = {v:k for k,v in subject2sample.items()}

filtered_sample_ids = [s for s in cidr_sample_ids if s in subject2sample]

for fh in glob.glob(f"{CIDR_GVCF_PREF}/*.gz"):
    for s in filtered_sample_ids:
        if subject2sample[s] in fh:
            smp2fh[s] = fh

SAMPLES = smp2fh.keys()

# BIN_SIZE = 5_000_000
# genome = pd.read_csv("hg38.genome", sep="\t")
# CHROMS = []
# STARTS, ENDS = [], []

# CHROM2INTERVAL = defaultdict(list)

# for i, row in genome.iterrows():
#     bins = np.arange(1, row["size"], BIN_SIZE)
#     starts, ends = bins[:-1], bins[1:]
#     if any([s in row["chrom"] for s in ("fix", "alt", "random", "Un")]): 
#         continue
#     if row["chrom"] != "chr21": continue
#     for start, end in zip(starts, ends):
#         if start == 1: continue
#         interval = f"{row['chrom']}_{start}_{end}"
#         CHROM2INTERVAL[row["chrom"]].append(interval)
#         # INTERVALS.append((row["chrom"], start, end))
#         # CHROMS.append(row["chrom"])
#         # STARTS.append(start)
#         # ENDS.append(end)

def get_orig_gvcf(wildcards):
    return smp2fh[wildcards.SAMPLE]

wildcard_constraints:
    INTERVAL = "chr[0-9]{1,2}_[0-9]+_[0-9]+",
    CHROM = "chr[0-9]{1,2}"

rule all:
    input:
        #"joint_called/vqsr/chr21.snp_indel.vcf.gz",
        "csv/slivar.chr21.tsv"

rule make_ped:
    input:
        orig_ped = f"{ORIG_PREF}/ped/2019-6-27-manual-checked-ceph.no_missing_samples.sorted.peddy_filtered.ped",
        cidr_ped = "Quinlan_Released_Data/Sample_Info/QuinlanNeklason_SIF.xlsx"
    output:
        ped = "joint.ped"
    run:
        import pandas as pd

        orig_ped = pd.read_csv(input.orig_ped, sep="\t", names=["FAMILY_ID", "SAMPLE_ID", "FATHER_ID", "MOTHER_ID", "SEX", "PHENOTYPE", "AGE", "DENOMINATOR"], dtype={"SAMPLE_ID": str})
        orig_ped = orig_ped[["FAMILY_ID", "SAMPLE_ID", "FATHER_ID", "MOTHER_ID", "SEX", "PHENOTYPE"]]
        orig_samples = orig_ped["SAMPLE_ID"].to_list()
        cidr_ped = pd.read_excel(input.cidr_ped, sheet_name="Sheet0", dtype={"Individual": str})
        cidr_ped = cidr_ped.rename(columns={"Family": "FAMILY_ID", "Individual": "SAMPLE_ID", "Father": "FATHER_ID", "Mother": "MOTHER_ID", "Sex": "SEX", "Phenotype": "PHENOTYPE"})
        cidr_ped = cidr_ped[["FAMILY_ID", "SAMPLE_ID", "FATHER_ID", "MOTHER_ID", "SEX"]]
        cidr_ped["PHENOTYPE"] = 0
        # ditch already existing samples
        cidr_ped = cidr_ped[~cidr_ped["SAMPLE_ID"].isin(orig_samples)]
        combined = pd.concat([orig_ped, cidr_ped])
        combined = combined[combined["FAMILY_ID"] == 1463]
        combined.to_csv(output.ped, index=False, sep="\t")


rule normalize_gvcf:
    input:
        gvcf = lambda wildcards: smp2fh[wildcards.SAMPLE]
    output:
        gvcf = "data/gvcf/{SAMPLE}.{CHROM}.normed.g.vcf.gz"
    threads: 4
    shell:
        """
        bcftools norm -r {wildcards.CHROM} \
                      --threads {threads} \
                      -m +any \
                      -Oz \
                      -o {output.gvcf} {input.gvcf}
        """


rule index_gvcf:
    input:
        gvcf = "data/gvcf/{SAMPLE}.{CHROM}.normed.g.vcf.gz"
    output:  "data/gvcf/{SAMPLE}.{CHROM}.normed.g.vcf.gz.tbi"
    shell:
        """
        module load bcftools
        
        bcftools index --tbi {input.gvcf}
        """


rule make_map:
    input:
        gvcfs = expand("data/gvcf/{SAMPLE}.{{CHROM}}.normed.g.vcf.gz", SAMPLE=SAMPLES),
        idxs = expand("data/gvcf/{SAMPLE}.{{CHROM}}.normed.g.vcf.gz.tbi", SAMPLE=SAMPLES)
    output: fh = "cohort.{CHROM}.sample_map"
    run:
        with open(output.fh, "w") as outfh:
            for fh in input.gvcfs:
                sample_id = fh.split("/")[-1].split(".")[0]
                print ("\t".join([sample_id, fh]), file=outfh)
        outfh.close()


rule add_interval_to_db:
    input:
        sample_map = "cohort.{CHROM}.sample_map"
    output:  
        db = directory("databases/{CHROM}_database/"),
    threads: 8
    resources:
        mem_mb = 32_000
    shell:
        """
        module load gatk/4.6

        echo {output.db}
        
        gatk --java-options '-Xmx24g -Xms24g' \
                GenomicsDBImport \
                --genomicsdb-workspace-path {output.db} \
                --batch-size 8 \
                -L {wildcards.CHROM} \
                --sample-name-map {input.sample_map} \
                --tmp-dir /scratch/ucgd/lustre-labs/quinlan/u1006375/gatk_tmp \
                --reader-threads 8

        """


rule genotype_gvcfs:
    input:
        database = "databases/{CHROM}_database/",
        ref = REF
    output:
        "joint_called/{CHROM}.vcf.gz"
    resources:
        mem_mb = 32_000
    shell:
        """
        module load gatk/4.6
        
        gatk --java-options '-Xmx24g -Xms24g' GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.database} \
            -O {output} \
            --intervals {wildcards.CHROM}
        """

rule make_site_level_vcf:
    input:
        "joint_called/{CHROM}.vcf.gz"
    output:
        "joint_called/sites/{CHROM}.sites.vcf.gz"
    resources:
        mem_mb = 32_000
    shell:
        """
        module load gatk/4.6

        gatk --java-options '-Xmx24g -Xms24g' MakeSitesOnlyVcf \
                -I {input} \
                -O {output}
        """

rule download_vqsr_indel_files:
    input:
    output:
        axiom = "vqsr_files/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
        mills = "vqsr_files/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        dbsnp = "vqsr_files/Homo_sapiens_assembly38.dbsnp138.vcf"
    shell:
        """
        wget -O {output.axiom} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
        wget -O {output.mills} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        wget -O {output.dbsnp} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf 
        """

rule compress_dbsnp:
    input:
        dbsnp = "vqsr_files/Homo_sapiens_assembly38.dbsnp138.vcf"
    output:
        dbsnp = "vqsr_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bgzip --threads {threads} {input.dbsnp}
        
        bcftools index --threads {threads} --tbi {output.dbsnp}
        """

rule download_vqsr_snp_files:
    input:
    output:
        thousand_genomes = "vqsr_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        hapmap = "vqsr_files/hapmap_3.3.hg38.vcf.gz",
        omni = "vqsr_files/1000G_omni2.5.hg38.vcf.gz",
    shell:
        """
        wget -O {output.thousand_genomes} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
        wget -O {output.hapmap} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
        wget -O {output.omni} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
        """

rule calculate_vqsr_indels:
    input:
        axiom = "vqsr_files/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
        mills = "vqsr_files/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        dbsnp = "vqsr_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        vcf = "joint_called/sites/{CHROM}.sites.vcf.gz"
    output:
        recal = "recal/{CHROM}.indel.recal",
        tranche = "tranche/{CHROM}.indel.tranches"
    resources: 
        mem_mb = 32_000
    shell:
        """

        module load gatk/4.6

        gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
             -V {input.vcf} \
             --trust-all-polymorphic \
             -tranche 100.0 \
             -tranche 99.95 \
             -tranche 99.9 \
             -tranche 99.5 \
             -tranche 99.0 \
             -tranche 97.0 \
             -tranche 96.0 \
             -tranche 95.0 \
             -tranche 94.0 \
             -tranche 93.5 \
             -tranche 93.0 \
             -tranche 92.0 \
             -tranche 91.0 \
             -tranche 90.0 \
             -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
             -mode INDEL \
             --max-gaussians 4 \
             -resource:mills,known=false,training=true,truth=true,prior=12 {input.mills} \
             -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {input.axiom} \
             -resource:dbsnp,known=true,training=false,truth=false,prior=2 {input.dbsnp} \
             -O {output.recal} \
             --tranches-file {output.tranche}
        """

rule calculate_vqsr_snps:
    input:
        dbsnp = "vqsr_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        thousand_genomes = "vqsr_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        hapmap = "vqsr_files/hapmap_3.3.hg38.vcf.gz",
        omni = "vqsr_files/1000G_omni2.5.hg38.vcf.gz",
        vcf = "joint_called/sites/{CHROM}.sites.vcf.gz"
    output:
        recal = "recal/{CHROM}.snp.recal",
        tranche = "tranche/{CHROM}.snp.tranches"
    shell:
        """
        module load gatk/4.6


        gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
                -V {input.vcf} \
                --trust-all-polymorphic \
                -tranche 100.0 \
                -tranche 99.95 \
                -tranche 99.9 \
                -tranche 99.8 \
                -tranche 99.6 \
                -tranche 99.5 \
                -tranche 99.4 \
                -tranche 99.3 \
                -tranche 99.0 \
                -tranche 98.0 \
                -tranche 97.0 \
                -tranche 90.0 \
                -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
                -mode SNP \
                --max-gaussians 6 \
                -resource:hapmap,known=false,training=true,truth=true,prior=15 {input.hapmap} \
                -resource:omni,known=false,training=true,truth=true,prior=12 {input.omni} \
                -resource:1000G,known=false,training=true,truth=false,prior=10 {input.thousand_genomes} \
                -resource:dbsnp,known=true,training=false,truth=false,prior=7 {input.dbsnp} \
                -O {output.recal} \
                --tranches-file {output.tranche}
        """

rule apply_vqsr_indels:
    input:
        vcf = "joint_called/{CHROM}.vcf.gz",
        recal = "recal/{CHROM}.indel.recal",
        tranche = "tranche/{CHROM}.indel.tranches"
    output:
        vcf = "joint_called/vqsr/{CHROM}.vcf.gz"
    resources:
        mem_mb = 8_000
    shell:
        """

        module load gatk/4.6

        gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
                -V {input.vcf} \
                --recal-file {input.recal} \
                --tranches-file {input.tranche} \
                --truth-sensitivity-filter-level 99.7 \
                --create-output-variant-index true \
                -mode INDEL \
                -O {output.vcf}
        """

rule apply_vqsr_snps:
    input:
        vcf = "joint_called/vqsr/{CHROM}.vcf.gz",
        recal = "recal/{CHROM}.snp.recal",
        tranche = "tranche/{CHROM}.snp.tranches"
    output:
        vcf = "joint_called/vqsr/{CHROM}.snp_indel.vcf.gz"
    resources:
        mem_mb = 8_000
    shell:
        """
        module load gatk/4.6

        gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --tranches-file {input.tranche} \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            -mode SNP \
            -O {output.vcf}
        """



rule merge_genotyped_vcfs:
    input:
        vcfs = lambda wildcards: expand("joint_called/vqsr/{CHROM}.snp_indel.vcf.gz", CHROM=CHROMS),
    output:
        vcf = "joint_called/merged.vqsr.vcf.gz"
    threads: 8
    resources:
        mem_mb = 32_000
    shell:
        """
        module load bcftools

        bcftools concat --threads {threads} -Oz -o {output.vcf} {input.vcfs}
        """

rule find_denovos:
    input:
        vcf = "joint_called/vqsr/{CHROM}.snp_indel.vcf.gz",
        ped = "joint.ped"
    output:
        vcf = "joint_called/slivar/{CHROM}.filtered.bcf"
    shell:
        """
        slivar expr \
            --pass-only \
            --vcf {input.vcf} \
            --ped {input.ped} \
            --out-vcf {output.vcf} \
            --info 'variant.call_rate > 0.9 && !variant.is_multiallelic' \
            --trio 'denovo:kid.het && mom.hom_ref && dad.hom_ref \
                            && kid.AB > 0.2 && kid.AB < 0.8 \
                            && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 \
                            && kid.DP >= 10 && mom.DP >= 10 && dad.DP >= 10' \
                          
        """

rule output_denovo_tsv:
    input:
        vcf = "joint_called/slivar/{CHROM}.filtered.bcf",
        ped = "joint.ped"
    output:
        tsv = "csv/slivar.{CHROM}.tsv"
    shell:
        """
        slivar tsv -p {input.ped} \
            -s denovo \
            {input.vcf} > {output.tsv}
        """