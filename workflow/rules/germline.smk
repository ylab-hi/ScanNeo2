rule picard_create_dict:
    input:
        "resources/refs/genome.fasta"
    output:
        "resources/refs/genome.dict"
    log:
        "logs/picard/create_dict.log"
    params:
        extra="",  # optional: extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/picard/createsequencedictionary"


# download training sets for calling high confidence variants
# see https://gatk.broadinstitute.org/hc/en-us/articles/4402736812443-Which-training-sets-arguments-should-I-use-for-running-VQSR-
rule gatk_vqsr_training_sets:
    output:
        snp_hapmap="resources/vqsr/hapmap_3.3.hg38.vcf.gz",
        snp_hapmap_idx="resources/vqsr/hapmap_3.3.hg38.vcf.gz.tbi",
        snp_omni="resources/vqsr/1000G_omni2.5.hg38.vcf.gz",
        snp_omni_idx="resources/vqsr/1000G_omni2.5.hg38.vcf.gz.tbi",
        snp_1000g="resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        snp_1000g_idx="resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
        snp_dbSNP="resources/vqsr/dbSNP_b150.vcf.gz",
        snp_dbSNP_idx="resources/vqsr/dbSNP_b150.vcf.gz.tbi",
        indel_mills="resources/vqsr/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        indel_mills_idx="resources/vqsr/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"

    shell:
        """
        curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz -o resources/vqsr/hapmap_3.3.hg38.vcf.gz
        curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi -o resources/vqsr/hapmap_3.3.hg38.vcf.gz.tbi
        curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz -o resources/vqsr/1000G_omni2.5.hg38.vcf.gz
        curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi -o resources/vqsr/1000G_omni2.5.hg38.vcf.gz.tbi
        curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz -o resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz
        curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi -o resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
        curl -L https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/00-common_all.vcf.gz  -o resources/vqsr/dbSNP_b150.vcf.gz
        curl -L https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/00-common_all.vcf.gz.tbi  -o resources/vqsr/dbSNP_b150.vcf.gz.tbi
        curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -o resources/vqsr/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        curl -L https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi -o resources/vqsr/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
        """

# do a first round of variant calling on original, unrecalibrated data
rule htc_first:
    input:
        # single or list of bam files
        bam="results/{sample}/rnaseq/align/realigned.bam",
        ref="resources/refs/genome.fasta",
        known="resources/vqsr/dbSNP_b150.vcf.gz"  # optional
    output:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.1rd.vcf"
    log:
        "logs/{sample}/gatk/haplotypecaller/1rd.log",
    params:
        extra="", 
        java_opts="",
    threads: config['threads']
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/haplotypecaller"


# recalibrate variants (SNP)
rule htc_first_snp_recal:
    input:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.1rd.vcf",
        ref="resources/refs/genome.fasta",
        hapmap="resources/vqsr/hapmap_3.3.hg38.vcf.gz",
        omni="resources/vqsr/1000G_omni2.5.hg38.vcf.gz",
        g1k="resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="resources/vqsr/dbSNP_b150.vcf.gz",
    output:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.1rd.snp.recal.vcf",
        idx="results/{sample}/rnaseq/indel/htcaller/variants.1rd.snp.recal.vcf.idx",
        tranches="results/{sample}/rnaseq/indel/htcaller/variants.1rd.snp.all.tranches"
    log:
        "logs/{sample}/gatk/vqsr/recal.first.snp"

    params:
        mode="SNP", 
        resources={
            "hapmap": {"known": False, "training": True, "truth": True, "prior": 15.0},
            "omni": {"known": False, "training": True, "truth": True, "prior": 12.0},
            "g1k": {"known": False, "training": True, "truth": True, "prior": 10.0},
            "dbsnp": {"known": False, "training": False, "truth": False, "prior": 2.0},
        },
        annotation=["MQ", "QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR", "DP"],
        extra=""
    threads: config['threads']
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/variantrecalibrator"


rule htc_first_snp_apply_vqsr:
    input:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.1rd.vcf",
        recal="results/{sample}/rnaseq/indel/htcaller/variants.1rd.snp.recal.vcf",
        tranches="results/{sample}/rnaseq/indel/htcaller/variants.1rd.snp.all.tranches",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.1rd.snp.filt.vcf"
    log:
        "logs/{sample}/gatk/vqsr/apply.1rd.snp.log"
    params:
        mode="SNP",  # set mode, must be either SNP, INDEL or BOTH
        extra="--truth-sensitivity-filter-level 99.5",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/applyvqsr"


# repeat Variant Quality Score Recalibration for indels
use rule htc_first_snp_recal as htc_first_indel_recal with:
    output:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.1rd.indel.recal.vcf",
        idx="results/{sample}/rnaseq/indel/htcaller/variants.1rd.indel.recal.vcf.idx",
        tranches="results/{sample}/rnaseq/indel/htcaller/variants.1rd.indel.all.tranches"
    log:
        "logs/{sample}/gatk/vqsr/recal_first_snp.log"
    params:
        mode="INDEL", 
        resources={
            "hapmap": {"known": False, "training": True, "truth": True, "prior": 15.0},
            "omni": {"known": False, "training": True, "truth": True, "prior": 12.0},
            "g1k": {"known": False, "training": True, "truth": True, "prior": 10.0},
            "dbsnp": {"known": False, "training": False, "truth": False, "prior": 2.0},
        },
        annotation=["MQ", "QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR", "DP"],
        extra=""


use rule htc_first_snp_apply_vqsr as htc_first_indel_apply_vsqr with:
    input:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.1rd.vcf",
        recal="results/{sample}/rnaseq/indel/htcaller/variants.1rd.indel.recal.vcf",
        tranches="results/{sample}/rnaseq/indel/htcaller/variants.1rd.indel.all.tranches",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.1rd.indel.filt.vcf"
    output:
    log:
        "logs/{sample}/gatk/vqsr/apply.final.indel.log"
    params:
        mode="INDEL", 
        extra="--truth-sensitivity-filter-level 99.0",  # optional


rule gatk_baserecalibrator:
    input:
        bam="results/{sample}/rnaseq/align/realigned.bam",
        ref="resources/refs/genome.fasta",
        dict="resources/refs/genome.dict",
        known=["results/{sample}/rnaseq/indel/htcaller/variants.1rd.indel.filt.vcf", 
            "results/{sample}/rnaseq/indel/htcaller/variants.1rd.snp.filt.vcf"]
    output:
        recal_table="results/{sample}/rnaseq/indel/htcaller/variants.1rd.baserecal.grp"
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=10240,
    wrapper:
        "v1.31.1/bio/gatk/baserecalibrator"


rule gatk_applybqsr:
    input:
        bam="results/{sample}/rnaseq/align/realigned.bam",
        ref="resources/refs/genome.fasta",
        dict="resources/refs/genome.dict",
        recal_table="results/{sample}/rnaseq/indel/htcaller/variants.1rd.baserecal.grp"
    output:
        bam="results/{sample}/rnaseq/indel/htcaller/variants.1rd.baserecal.bam"
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
        embed_ref=True,  # embed the reference in cram output
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/applybqsr"

rule htcaller_main:
    input:
        # single or list of bam files
        bam="results/{sample}/rnaseq/indel/htcaller/variants.1rd.baserecal.bam",
        ref="resources/refs/genome.fasta",
        known="resources/vqsr/dbSNP_b150.vcf.gz"  # optional
    output:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.final.vcf"
    log:
        "logs/gatk/haplotypecaller/{sample}.1rd.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    threads: config['threads']
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/haplotypecaller"


use rule htc_first_snp_recal as htc_final_snp_recal with:
    input:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.final.vcf",
        ref="resources/refs/genome.fasta",
        hapmap="resources/vqsr/hapmap_3.3.hg38.vcf.gz",
        omni="resources/vqsr/1000G_omni2.5.hg38.vcf.gz",
        g1k="resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="resources/vqsr/dbSNP_b150.vcf.gz",
    output:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.final.snp.recal.vcf",
        idx="results/{sample}/rnaseq/indel/htcaller/variants.final.snp.recal.vcf.idx",
        tranches="results/{sample}/rnaseq/indel/htcaller/variants.final.snp.all.tranches"
    log:
        "logs/{sample}/gatk/vqsr/recal.final.snp.log"


use rule htc_first_snp_apply_vqsr as htc_final_snp_apply_vsqr with:
    input:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.final.vcf",
        recal="results/{sample}/rnaseq/indel/htcaller/variants.final.snp.recal.vcf",
        tranches="results/{sample}/rnaseq/indel/htcaller/variants.final.snp.all.tranches",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/variants/germ.snvs.vcf"
    log:
        "logs/{sample}/gatk/vqsr/apply.final.snp.log"


use rule htc_first_snp_recal as htc_final_indel_recal with:
    input:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.final.vcf",
        ref="resources/refs/genome.fasta",
        hapmap="resources/vqsr/hapmap_3.3.hg38.vcf.gz",
        omni="resources/vqsr/1000G_omni2.5.hg38.vcf.gz",
        g1k="resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="resources/vqsr/dbSNP_b150.vcf.gz",
    output:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.final.indel.recal.vcf",
        idx="results/{sample}/rnaseq/indel/htcaller/variants.final.indel.recal.vcf.idx",
        tranches="results/{sample}/rnaseq/indel/htcaller/variants.final.indel.all.tranches"
    log:
        "logs/{sample}/gatk/vqsr/recal.final.indel.log"
    params:
        mode="INDEL",  
        resources={
            "hapmap": {"known": False, "training": True, "truth": True, "prior": 15.0},
            "omni": {"known": False, "training": True, "truth": True, "prior": 12.0},
            "g1k": {"known": False, "training": True, "truth": True, "prior": 10.0},
            "dbsnp": {"known": False, "training": False, "truth": False, "prior": 2.0},
        },
        annotation=["MQ", "QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR", "DP"],
        extra=""
        

use rule htc_first_snp_apply_vqsr as htc_final_indel_apply_vsqr with:
    input:
        vcf="results/{sample}/rnaseq/indel/htcaller/variants.final.vcf",
        recal="results/{sample}/rnaseq/indel/htcaller/variants.final.indel.recal.vcf",
        tranches="results/{sample}/rnaseq/indel/htcaller/variants.final.indel.all.tranches",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/variants/germ.indel.vcf"
    log:
        "logs/{sample}/gatk/vqsr/apply.final.indel.log"
    params:
        mode="INDEL", 
        extra="--truth-sensitivity-filter-level 99.0",  # optional
