import pandas as pd
import glob
from collections import defaultdict
import numpy as np
import csv
import json
from collections import namedtuple

### USAGE NOTE: run with -R add_interval_to_db if APPEND is true

# constrain wildcards
wildcard_constraints:
    INTERVAL = "chr[0-9]{1,2}_[0-9]+_[0-9]+",
    CHROM = "chr[0-9]{1,2}|chrX|chrY",
    START = "[0-9]+",
    END = "[0-9]+"

# read in the samples we want to use for this run
sample_fh = config["sample_path"]
SAMPLES = []
with open(sample_fh, "r") as infh:
    for l in infh:
        SAMPLES.append(l.strip())


# figure out whether we want to append to genomics-db-workspaces or
# instantiate from scratch
APPEND = bool(config["append"])

EXISTING_SAMPLES = []
# figure out if any of these samples are already in the genomics-db
for json_fh in glob.glob("genomics_db/**/callset.json"):
    with open(json_fh) as f:
        dicts = json.load(f)
        for d in dicts["callsets"]:
            EXISTING_SAMPLES.append(d["sample_name"])

OVERLAPPING_SAMPLES = [s for s in SAMPLES if s in EXISTING_SAMPLES]
if any([s in EXISTING_SAMPLES for s in SAMPLES]):
    print (f"### WARNING!!! {len(OVERLAPPING_SAMPLES)} of the specified samples have already been joint-genotyped! ###\n{','.join(OVERLAPPING_SAMPLES)}")

# map samples to GVCFs
SMP2GVCF = {}
for json_fh in glob.glob("json/gvcf_mapping.json"):
    with open(json_fh) as f:
        dicts = json.load(f)
        for d in dicts:
            sample = d["sample"]
            gvcf_fh = d["gvcf_fh"]
            SMP2GVCF[sample] = gvcf_fh

REF = "/scratch/ucgd/lustre/common/data/Reference/homo_sapiens/GRCh38/primary_assembly_decoy_phix_masked.fa"

db_interval = namedtuple("interval", ["chrom", "start", "end"])
INTERVALS = []
with open("gatk_intervals/wgs_calling_regions.hg38.interval_list.1Mbp.bed", "rt") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    for l in csvf:
        chrom, start, end = l
        if chrom != "chr22": continue
        INTERVALS.append(db_interval(chrom, start, end))

rule all:
    input:
        expand("genomics_db/{CHROM}/{START}_{END}_database/completed.txt", zip, CHROM=[nt.chrom for nt in INTERVALS], START=[nt.start for nt in INTERVALS], END=[nt.end for nt in INTERVALS])
        # "genomics_db/chr22/16154319_16279672_database/completed.txt"
        
rule normalize_gvcf:
    input:
        gvcf = lambda wildcards: SMP2GVCF[wildcards.SAMPLE],
    output:
        gvcf = temp("data/gvcf/{SAMPLE}.{CHROM}.normed.g.vcf.gz")
    resources:
        mem_mb = 32_000
    threads: 4
    shell:
        """
        module load bcftools

        bcftools norm --regions {wildcards.CHROM} \
                      --threads {threads} \
                      -m +any \
                      -Oz \
                      -o {output.gvcf} \
                      {input.gvcf}
        """


rule index_normalized_gvcf:
    input:
        gvcf = "data/gvcf/{SAMPLE}.{CHROM}.normed.g.vcf.gz"
    output:
        idx = temp("data/gvcf/{SAMPLE}.{CHROM}.normed.g.vcf.gz.tbi")
    shell:
        """
        module load gatk/4.6

        gatk IndexFeatureFile -I {input.gvcf}
        """


rule remove_mnps:
    input:
        gvcf = "data/gvcf/{SAMPLE}.{CHROM}.normed.g.vcf.gz",
        gvcf_idx = "data/gvcf/{SAMPLE}.{CHROM}.normed.g.vcf.gz.tbi",
        ref = REF,
    output:
        gvcf = "data/gvcf/{SAMPLE}.{CHROM}.normed.no_mnps.g.vcf.gz"
    resources:
        mem_mb = 32_000
    shell:
        """
        module load gatk/4.6

        gatk --java-options '-Xmx24g -Xms24g' \
                SelectVariants \
                -V {input.gvcf} \
                -R {input.ref} \
                --select-type-to-exclude MNP \
                -O {output.gvcf} \
                --ignore-non-ref-in-types
        """


rule index_gvcf:
    input:
        gvcf = "data/gvcf/{SAMPLE}.{CHROM}.normed.no_mnps.g.vcf.gz"
    output: "data/gvcf/{SAMPLE}.{CHROM}.normed.no_mnps.g.vcf.gz.tbi"

    shell:
        """
        module load gatk/4.6
        
        gatk IndexFeatureFile -I {input.gvcf}
        """


rule make_map:
    input:
        gvcfs = expand("data/gvcf/{SAMPLE}.{{CHROM}}.normed.no_mnps.g.vcf.gz", SAMPLE=SAMPLES),
        idxs = expand("data/gvcf/{SAMPLE}.{{CHROM}}.normed.no_mnps.g.vcf.gz.tbi", SAMPLE=SAMPLES)
    output: fh = "sample_maps/{CHROM}.sample_map"
    run:
        with open(output.fh, "w") as outfh:
            for fh in input.gvcfs:
                sample_id = fh.split("/")[-1].split(".")[0]
                print ("\t".join([sample_id, fh]), file=outfh)
        outfh.close()


rule add_interval_to_db:
    input:
        sample_map = "sample_maps/{CHROM}.sample_map",
        gvcfs = expand("data/gvcf/{SAMPLE}.{{CHROM}}.normed.no_mnps.g.vcf.gz", SAMPLE=SAMPLES),
        idxs = expand("data/gvcf/{SAMPLE}.{{CHROM}}.normed.no_mnps.g.vcf.gz.tbi", SAMPLE=SAMPLES)
    output:  
        # snakemake will remove this file before running the rule, but that's ok since the file is meaningless
        completed = "genomics_db/{CHROM}/{START}_{END}_database/completed.txt"
    threads: 16
    params:
        interval_arg = lambda wildcards: "" if APPEND else f"-L {wildcards.CHROM}:{wildcards.START}-{wildcards.END}",
        workspace_arg = "--genomicsdb-update-workspace-path" if APPEND else "--genomicsdb-workspace-path",
        db_dir = "genomics_db/{CHROM}/{START}_{END}_database/"
    resources:
        mem_mb = 32_000
    script:
        "scripts/add_interval_to_db.sh"


