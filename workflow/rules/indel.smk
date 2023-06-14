import os
from snakemake.remote import HTTP

rule realign:
    input:
        idx = multiext("resources/refs/bwa/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        bam = "results/{sample}/rnaseq/preproc/post/aln.bam"
    output:
        "results/{sample}/rnaseq/preproc/post/realn.bam"
    log:
        "logs/bwa_realign/{sample}.log"
    conda:
        "../envs/realign.yml",
    threads: config['threads']
    shell:
        """
        samtools collate -Oun128 {input.bam} \
        | samtools fastq -OT RG,BC - \
        | bwa mem -pt8 -CH <(samtools view -H {input.bam}|grep ^@RG) resources/refs/bwa/genome - \
        | samtools sort -@{threads} -m4g - \
        | samtools addreplacerg -r '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}' -o {output} -
        """

rule realign_index:
    input:
        "results/{sample}/rnaseq/preproc/post/realn.bam"
    output:
        "results/{sample}/rnaseq/preproc/post/realn.bam.bai",
    log:
        "logs/samtools_index_{sample}.log"
    params:
        extra="",  # optional params string
    threads: config['threads']  # This value - 1 will be sent to -@
    wrapper:
        "v1.26.0/bio/samtools/index"
        

rule transindel_build:
    input:
        bam = "results/{sample}/rnaseq/preproc/post/realn.bam",
        idx = "results/{sample}/rnaseq/preproc/post/realn.bam.bai"
    output:
        "results/{sample}/rnaseq/indel/transindel/build.bam"
    log:
        "logs/transindel/build/{sample}.log"
    conda:
        "../envs/transindel.yml"
    shell:
        """
        python3 workflow/scripts/transIndel/transIndel_build_RNA.py \
        -i {input.bam} \
        -o {output} \
        -r resources/refs/genome.fasta \
        -g resources/refs/genome.gtf > {log} 
        """

# builds index file for 
rule transindel_build_index:
    input:
        "results/{sample}/rnaseq/indel/transindel/build.bam"
    output:
        "results/{sample}/rnaseq/indel/transindel/build.bam.bai"
    log:
        "logs/samtools/index/transindel_build/{sample}.log"
    conda:
        "../envs/transindel.yml"
    shell:
        "samtools index {input} -o {output} > {log}"
    

rule transindel_call:
    input:
        bam = "results/{sample}/rnaseq/indel/transindel/build.bam",
        bai = "results/{sample}/rnaseq/indel/transindel/build.bam.bai"
    output:
        "results/{sample}/rnaseq/indel/transindel/call.indel.vcf"
    log:
        "logs/transindel/call/{sample}.log"
    conda:
        "../envs/transindel.yml"
    params:
        mapq=config['mapq']
    shell:
        """
        python workflow/scripts/transIndel/transIndel_call.py \
        -i {input.bam} \
        -l 10 \
        -o results/{wildcards.sample}/rnaseq/indel/transindel/call.indel.vcf \
        -m {params}
        """

# resove alleles and remove PCR slippage
rule slippage_removal:
    input:
        "results/{sample}/rnaseq/indel/transindel/call.indel.vcf"
    output:
        "results/{sample}/rnaseq/indel/transindel.vcf"
    log:
        "logs/indel/sliprem{sample}.log"
    conda:
        "../envs/transindel.yml"
    shell:
        """
            python3 workflow/scripts/slippage_removal.py \
            resources/refs/genome.fasta {input} {output} > {log} 
        """

rule gatk_mutect2:
    input:
        fasta="resources/refs/genome.fasta",
        map="results/indel/realign/{sample}.bam",
    output:
        vcf="results/indel/mutect2/{sample}.vcf",
        bam="results/indel/mutect2/{sample}.bam",
    message:
        "Mutect2 with {wildcards.sample}"
    threads: config['threads']
    resources:
        mem_mb=1024,
    params:
        extra="--min-base-quality-score 20",
    log:
        "logs/mutect_{sample}.log",
    wrapper:
        "v1.31.1/bio/gatk/mutect"


rule gatk_filtermutectcalls:
    input:
        vcf="results/indel/mutect2/{sample}.vcf",
        ref="resources/refs/genome.fasta",
        bam="results/indel/mutect2/{sample}.bam",
        # intervals="intervals.bed",
        # contamination="", # from gatk CalculateContamination
        # segmentation="", # from gatk CalculateContamination
        # f1r2="", # from gatk LearnReadOrientationBias
    output:
        vcf="results/indel/mutect2/{sample}_filtered.vcf"
    log:
        "logs/gatk/filter/snvs_{sample}.log",
    params:
        extra="--max-alt-allele-count 3",
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/filtermutectcalls"


rule gatk_select_SNPs:
    input:
        vcf="results/indel/mutect2/{sample}_filtered.vcf",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/indel/mutect2/{sample}.snvs.vcf",
    log:
        "logs/gatk/select/{sample}.snvs.log",
    params:
        extra="--select-type-to-include SNP",  # optional filter arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/selectvariants"


rule gatk_select_Indels:
    input:
        vcf="results/indel/mutect2/{sample}_filtered.vcf",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/indel/mutect2/{sample}.indels.vcf",
    log:
        "logs/gatk/select/{sample}.snvs.log",
    params:
        extra="--select-type-to-include INDEL",  # optional filter arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/selectvariants"


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
rule haplotype_caller_first_round:
    input:
        # single or list of bam files
        bam="results/indel/realign/{sample}.bam",
        ref="resources/refs/genome.fasta",
        known="resources/vqsr/dbSNP_b150.vcf.gz"  # optional
    output:
        vcf="results/indel/haplotypecaller/{sample}/{sample}.1rd.vcf"
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


# recalibrate variants (SNP)
rule gatk_variant_recal_snp:
    input:
        vcf="results/indel/haplotypecaller/{sample}/{sample}.1rd.vcf",
        ref="resources/refs/genome.fasta",
        hapmap="resources/vqsr/hapmap_3.3.hg38.vcf.gz",
        omni="resources/vqsr/1000G_omni2.5.hg38.vcf.gz",
        g1k="resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="resources/vqsr/dbSNP_b150.vcf.gz",
        
    output:
        vcf="results/indel/haplotypecaller/{sample}/{sample}.1rd.snp.recal.vcf",
        idx="results/indel/haplotypecaller/{sample}/{sample}.1rd.snp.recal.vcf.idx",
        tranches="results/indel/haplotypecaller/{sample}/{sample}.1rd.snp.all.tranches"

    log:
        "logs/gatk/variantrecalibrator/{sample}.log",

    params:
        mode="SNP",  # set mode, must be either SNP, INDEL or BOTH
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


rule gatk_apply_vqsr_snp:
    input:
        vcf="results/indel/haplotypecaller/{sample}/{sample}.1rd.vcf",
        recal="results/indel/haplotypecaller/{sample}/{sample}.1rd.snp.recal.vcf",
        tranches="results/indel/haplotypecaller/{sample}/{sample}.1rd.snp.all.tranches",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/indel/haplotypecaller/{sample}/{sample}.1rd.snp.filt.vcf",
    log:
        "logs/gatk/applyvqsr/{sample}_snp.log",
    params:
        mode="SNP",  # set mode, must be either SNP, INDEL or BOTH
        extra="--truth-sensitivity-filter-level 99.5",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "master/bio/gatk/applyvqsr"


# recalibrate variants (indel)
rule gatk_variant_recal_indel:
    input:
        vcf="results/indel/haplotypecaller/{sample}/{sample}.1rd.vcf",
        ref="resources/refs/genome.fasta",
        mills="resources/vqsr/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        dbsnp="resources/vqsr/dbSNP_b150.vcf.gz",
        
    output:
        vcf="results/indel/haplotypecaller/{sample}/{sample}.1rd.indel.recal.vcf",
        idx="results/indel/haplotypecaller/{sample}/{sample}.1rd.indel.recal.vcf.idx",
        tranches="results/indel/haplotypecaller/{sample}/{sample}.1rd.indel.all.tranches"

    log:
        "logs/gatk/variantrecalibrator/{sample}.log",

    params:
        mode="INDEL",  # set mode, must be either SNP, INDEL or BOTH
        resources={
            "mills": {"known": False, "training": True, "truth": True, "prior": 12.0},
            "dbsnp": {"known": True, "training": False, "truth": False, "prior": 2.0},
        },
        annotation=["MQ", "QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR", "DP"],
        extra="--max-gaussians 4"
    threads: config['threads']
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/variantrecalibrator"


rule gatk_apply_vqsr_indel:
    input:
        vcf="results/indel/haplotypecaller/{sample}/{sample}.1rd.vcf",
        recal="results/indel/haplotypecaller/{sample}/{sample}.1rd.indel.recal.vcf",
        tranches="results/indel/haplotypecaller/{sample}/{sample}.1rd.indel.all.tranches",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/indel/haplotypecaller/{sample}/{sample}.1rd.indel.filt.vcf",
    log:
        "logs/gatk/applyvqsr/{sample}_indel.log",
    params:
        mode="INDEL",  # set mode, must be either SNP, INDEL or BOTH
        extra="--truth-sensitivity-filter-level 99.0",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/applyvqsr"


rule gatk_baserecalibrator:
    input:
        bam="results/indel/realign/{sample}.bam",
        ref="resources/refs/genome.fasta",
        dict="resources/refs/genome.dict",
        known=["results/indel/haplotypecaller/{sample}/{sample}.1rd.indel.filt.vcf", "results/indel/haplotypecaller/{sample}/{sample}.1rd.snp.filt.vcf"]
    output:
        recal_table="results/indel/recal/{sample}.grp",
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/baserecalibrator"

#rule gatk_applybqsr:
#    input:
#        bam="results/indel/realign/{sample}.bam",
#        ref="resources/refs/genome.fasta",
#        dict="resources/refs/genome.dict",
#        recal_table="recal/{sample}.grp",
#    output:
#        bam="recal/{sample}.bam",
#    log:
#        "logs/gatk/gatk_applybqsr/{sample}.log",
#    params:
#        extra="",  # optional
#        java_opts="",  # optional
#        embed_ref=True,  # embed the reference in cram output
#    resources:
#        mem_mb=1024,
#    wrapper:
#        "v1.31.1/bio/gatk/applybqsr"



#rule picard_split_vcfs:
#    input:
#        "results/indel/haplotypecaller/{sample}.vcf"
#    output:
#        snp="results/indel/{sample}.snp.vcf",
#        indel="results/indel/{sample}.indel.vcf"
#    log:
#        "logs/spltvcfs{sample}.log"
#    conda:
#        "../envs/picard.yml"
#    shell:
#        """
#            picard SplitVcfs I={input} \
#            SNP_OUTPUT={output.snp} \
#            INDEL_OUTPUT={output.indel}
#            STRICT=false
#        """





