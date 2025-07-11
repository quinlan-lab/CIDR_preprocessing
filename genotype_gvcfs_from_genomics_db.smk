REF = "/scratch/ucgd/lustre/common/data/Reference/homo_sapiens/GRCh38/primary_assembly_decoy_phix_masked.fa"

CHROMS = list(map(str, range(1, 23)))
CHROMS = [f"chr{c}" for c in CHROMS] + ["chrX", "chrY"]

CHROM2INTERVAL = defaultdict(list)
with open("gatk_inteverals/wgs_calling_regions.hg38.interval_list.1Mbp.bed", "rt") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    for l in csvf:
        chrom, start, end = l
        CHROM2INTERVAL[chrom].append((start, end))


rule all:
    input: 
        expand("vcf/joint_called/{CHROM}.merged.vcf.gz", CHROM=CHROMS)


rule genotype_gvcfs:
    input:
        database = "databases/{CHROM}/{START}_{END}_database/",
        ref = REF
    output: "vcf/joint_called/{CHROM}/{START}.{END}.vcf.gz"
    resources:
        mem_mb = 32_000
    shell:
        """
        module load gatk/4.6
        
        gatk --java-options '-Xmx24g -Xms24g' GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.database} \
            -O {output} \
            --intervals {wildcards.CHROM}:{wildcards.START}-{wildcards.END}
        """


rule merge_genotyped_vcfs:
    input:
        vcfs = lambda wildcards: expand("vcf/joint_called/{{CHROM}}/{START}.{END}.vcf.gz", zip, START=[s for s,e in CHROM2INTERVAL[wildcards.CHROM]], END=[e for s,e in CHROM2INTERVAL[wildcards.CHROM]]),
    output:
        vcf = "vcf/joint_called/{CHROM}.merged.vcf.gz"
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
        vcf = "vcf/joint_called/{CHROM}.merged.vcf.gz",
        ped = "ped/joint.cidr.ceph.ped"
    output:
        vcf = "vcf/slivar/{CHROM}.filtered.bcf"
    shell:
        """
        slivar expr \
            --pass-only \
            --vcf {input.vcf} \
            --ped {input.ped} \
            --out-vcf {output.vcf} \
            --info '!variant.is_multiallelic' \
            --trio 'denovo:kid.het && mom.hom_ref && dad.hom_ref \
                            && kid.AB >= 0.2 && kid.AB <= 0.8 \
                            && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 \
                            && kid.DP >= 10 && mom.DP >= 10 && dad.DP >= 10' \
                          
        """

rule output_denovo_tsv:
    input:
        vcf = "vcf/slivar/{CHROM}.filtered.bcf",
        ped = "ped/joint.cidr.ceph.ped"
    output:
        tsv = "csv/slivar.{CHROM}.tsv"
    shell:
        """
        slivar tsv -p {input.ped} \
            -s denovo \
            {input.vcf} > {output.tsv}
        """