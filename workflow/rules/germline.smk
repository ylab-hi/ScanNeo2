# download training sets for calling high confidence variants
# see https://gatk.broadinstitute.org/hc/en-us/articles/4402736812443-Which-training-sets-arguments-should-I-use-for-running-VQSR-
rule get_gatk_vqsr_training_sets:
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
  message:
    "Downloading training sets for calling high confidence variants"
  log:
    "logs/gatk_get_training)_sets.log"
  conda:
    "../envs/basic.yml"
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
rule detect_variants_htc_first_round:
  input:
    # single or list of bam files
    bam="results/{sample}/{seqtype}/align/{group}_final_BWA.bam",
    ref="resources/refs/genome.fasta",
    known="resources/vqsr/dbSNP_b150.vcf.gz"  # optional
  output:
    vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.vcf"
  message:
    "First round of variant calling on original, unrecalibrated data on sample:{wildcards.sample} with group:{wildcards.group}" 
  log:
    "logs/{sample}/gatk/haplotypecaller/{seqtype}_{group}_1rd.log",
  params:
    extra="", 
    java_opts="",
  threads: config['threads']
  resources:
    mem_mb=1024,
  wrapper:
    "v1.31.1/bio/gatk/haplotypecaller"


# recalibrate variants (SNP)
rule recalibrate_variants_first_round:
  input:
    vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.vcf",
    ref="resources/refs/genome.fasta",
    hapmap="resources/vqsr/hapmap_3.3.hg38.vcf.gz",
    omni="resources/vqsr/1000G_omni2.5.hg38.vcf.gz",
    g1k="resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
    dbsnp="resources/vqsr/dbSNP_b150.vcf.gz",
  output:
    vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.recal.vcf",
    idx="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.recal.vcf.idx",
    tranches="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.tranches"
  message:
    "Recalibrate variants (SNP) on sample:{wildcards.sample} with group:{wildcards.group}"
  log:
    "logs/{sample}/gatk/vqsr/{seqtype}_{group}_variants.1rd.recal.log"
  params:
    mode="BOTH", 
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

rule apply_VQSR_SNVs_first_round:
  input:
      vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.vcf",
      recal="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.recal.vcf",
      tranches="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.tranches",
      ref="resources/refs/genome.fasta",
      dict="resources/refs/genome.dict"
  output:
      vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_snvs.1rd.flt.vcf"
  message:
    "Apply VQSR (SNP) on sample:{wildcards.sample} with group:{wildcards.group}"
  log:
      "logs/{sample}/gatk/vqsr/{seqtype}_{group}_apply.1rd.snvs.log"
  params:
      mode="SNP",  # set mode, must be either SNP, INDEL or BOTH
      extra="--truth-sensitivity-filter-level 99.5",  # optional
  resources:
      mem_mb=1024,
  wrapper:
      "v1.31.1/bio/gatk/applyvqsr"

rule apply_VSQR_INDEL_first_round:
  input:
    vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.vcf",
    recal="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.recal.vcf",
    tranches="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.tranches",
    ref="resources/refs/genome.fasta",
    dict="resources/refs/genome.dict"
  output:
        vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_indels.1rd.flt.vcf"
  message:
    "Apply VQSR (INDEL) on sample:{wildcards.sample} with group:{wildcards.group}"
  log:
      "logs/{sample}/gatk/vqsr/{seqtype}_{group}_apply.1rd.indels.log"
  params:
      mode="INDEL",  # set mode, must be either SNP, INDEL or BOTH
      extra="--truth-sensitivity-filter-level 99.5",  # optional
  resources:
      mem_mb=1024,
  wrapper:
      "v1.31.1/bio/gatk/applyvqsr"
    
rule recalibrate_bases:
  input:
    bam="results/{sample}/{seqtype}/align/{group}_final_BWA.bam",
    ref="resources/refs/genome.fasta",
    dict="resources/refs/genome.dict",
    known=["results/{sample}/{seqtype}/indel/htcaller/{group}_indels.1rd.flt.vcf", 
      "results/{sample}/{seqtype}/indel/htcaller/{group}_snvs.1rd.flt.vcf"]
  output:
    recal_table="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.baserecal.grp"
  log:
    "logs/{sample}/gatk/baserecalibrator/{seqtype}_{group}.log"
  params:
    extra="",  # optional
    java_opts="",  # optional
  resources:
    mem_mb=10240,
  wrapper:
    "v1.31.1/bio/gatk/baserecalibrator"

rule apply_BQSR:
  input:
    bam="results/{sample}/{seqtype}/align/{group}_final_BWA.bam",
    ref="resources/refs/genome.fasta",
    dict="resources/refs/genome.dict",
    recal_table="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.baserecal.grp"
  output:
    bam="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.baserecal.bam"
  log:
    "logs/{sample}/gatk/gatk_applybqsr/{seqtype}_{group}.log",
  params:
    extra="",  # optional
    java_opts="",  # optional
    embed_ref=True,  # embed the reference in cram output
  resources:
    mem_mb=1024,
  wrapper:
    "v1.31.1/bio/gatk/applybqsr"

rule detect_variants_htc_final_round:
    input:
      # single or list of bam files
      bam="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.baserecal.bam",
      ref="resources/refs/genome.fasta",
      known="resources/vqsr/dbSNP_b150.vcf.gz"  # optional
    output:
      vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.vcf"
    log:
      "logs/{sample}/haplotypecaller/{seqtype}_{group}.1rd.log",
    params:
      extra="",  # optional
      java_opts="",  # optional
    threads: config['threads']
    resources:
      mem_mb=1024,
    wrapper:
      "v1.31.1/bio/gatk/haplotypecaller"

# recalibrate variants 
rule recalibrate_variants_final_round:
  input:
    vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.vcf",
    ref="resources/refs/genome.fasta",
    hapmap="resources/vqsr/hapmap_3.3.hg38.vcf.gz",
    omni="resources/vqsr/1000G_omni2.5.hg38.vcf.gz",
    g1k="resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
    dbsnp="resources/vqsr/dbSNP_b150.vcf.gz",
  output:
    vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.recal.vcf",
    idx="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.recal.vcf.idx",
    tranches="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.tranches"
  message:
    "Recalibrate variants (SNP) on sample:{wildcards.sample} with group:{wildcards.group}"
  log:
    "logs/{sample}/gatk/vqsr/{seqtype}_{group}_final.variants.recal.log"

  params:
    mode="BOTH", 
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


rule apply_VQSR_SNVs_final:
  input:
    vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.vcf",
    recal="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.recal.vcf",
    tranches="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.tranches",
    ref="resources/refs/genome.fasta",
    dict="resources/refs/genome.dict"
  output:
    vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_snvs.final.flt.vcf"
  message:
    "Apply VQSR (SNP) on sample:{wildcards.sample} with group:{wildcards.group}"
  log:
    "logs/{sample}/gatk/vqsr/{seqtype}_{group}_apply.final.snvs.log"
  params:
    mode="BOTH",  # set mode, must be either SNP, INDEL or BOTH
    extra="--truth-sensitivity-filter-level 99.5",  # optional
  resources:
      mem_mb=1024,
  wrapper:
      "v1.31.1/bio/gatk/applyvqsr"

rule apply_VSQR_INDELS_final:
  input:
    vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.vcf",
    recal="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.recal.vcf",
    tranches="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.final.tranches",
    ref="resources/refs/genome.fasta",
  output:
    vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_indel.final.flt.vcf"
  message:
    "Apply VQSR (INDEL) on sample:{wildcards.sample} with group:{wildcards.group}"
  log:
    "logs/{sample}/gatk/vqsr/{seqtype}_{group}_apply.final.indel.log"
  params:
    mode="INDEL",  # set mode, must be either SNP, INDEL or BOTH
    extra="--truth-sensitivity-filter-level 99.5",  # optional
  resources:
    mem_mb=1024,
  wrapper:
    "v1.31.1/bio/gatk/applyvqsr"

