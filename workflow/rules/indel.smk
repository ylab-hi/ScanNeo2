import os
from snakemake.remote import HTTP

rule transindel_build:
    input:
        bam = "results/{sample}/rnaseq/align/realigned.bam",
        idx = "results/{sample}/rnaseq/align/realigned.bam.bai"
    output:
        bam="results/{sample}/rnaseq/indel/transindel/build.bam",
        idx="results/{sample}/rnaseq/indel/transindel/build.bam.bai"
    log:
        "logs/transindel/build/{sample}.log"
    conda:
        "../envs/transindel.yml"
    shell:
        """
        python3 workflow/scripts/transIndel/transIndel_build_RNA.py \
        -i {input.bam} \
        -o {output.bam} \
        -r resources/refs/genome.fasta \
        -g resources/refs/genome.gtf > {log} 
        samtools index {output.bam} -o {output.idx} >> {log} 
        """

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
        -o results/{wildcards.sample}/rnaseq/indel/transindel/call \
        -m {params}
        """

# resove alleles and remove PCR slippage
rule slippage_removal:
    input:
        "results/{sample}/rnaseq/indel/transindel/call.indel.vcf"
    output:
        "results/{sample}/variants/long.indel.vcf"
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
        map="results/{sample}/rnaseq/align/realigned.bam"
    output:
        vcf="results/{sample}/rnaseq/indel/mutect2/variants.vcf",
        bam="results/{sample}/rnaseq/indel/mutect2/variants.bam",
    message:
        "Detection of somatic SNVs/Indels with Mutect2 on sample {wildcards.sample}"
    threads: config['threads']
    resources:
        mem_mb=10024,
    params:
        extra="",
    log:
        "logs/gatk/mutect2/{sample}.log",
    wrapper:
        "v1.31.1/bio/gatk/mutect"


rule gatk_filtermutectcalls:
    input:
        vcf="results/{sample}/rnaseq/indel/mutect2/variants.vcf",
        bam="results/{sample}/rnaseq/indel/mutect2/variants.bam",
        ref="resources/refs/genome.fasta",
        # intervals="intervals.bed",
        # contamination="", # from gatk CalculateContamination
        # segmentation="", # from gatk CalculateContamination
        # f1r2="", # from gatk LearnReadOrientationBias
    output:
        vcf="results/{sample}/rnaseq/indel/mutect2/variants_flt.vcf"
    log:
        "logs/gatk/filtermutect/{sample}.log",
    params:
        extra="--max-alt-allele-count 3",
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/filtermutectcalls"


rule gatk_select_SNPs:
    input:
        vcf="results/{sample}/rnaseq/indel/mutect2/variants_flt.vcf",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/variants/snvs.vcf"
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
        vcf="results/{sample}/rnaseq/indel/mutect2/variants_flt.vcf",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/variants/short.indel.vcf"
    log:
        "logs/gatk/select/{sample}.indel.log",
    params:
        extra="--select-type-to-include INDEL",  # optional filter arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/selectvariants"




# recalibrate variants (SNP)
#rule gatk_variant_recal_snp2:
#    input:
#        vcf="results/{sample}/rnaseq/indel/htcaller/variants.1rd.vcf",
#        ref="resources/refs/genome.fasta",
#        hapmap="resources/vqsr/hapmap_3.3.hg38.vcf.gz",
#        omni="resources/vqsr/1000G_omni2.5.hg38.vcf.gz",
#        g1k="resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
#        dbsnp="resources/vqsr/dbSNP_b150.vcf.gz",
        
#    output:
#        vcf="results/{sample}/rnaseq/indel/htcaller/variants.1rd.snp.recal.vcf",
#        idx="results/{sample}/rnaseq/indel/htcaller/variants.1rd.snp.recal.vcf.idx",
#        tranches="results/{sample}/rnaseq/indel/htcaller/variants.1rd.snp.all.tranches"
#    log:
#        "logs/gatk/variantrecalibrator/{sample}.log",

#    params:
#        mode="SNP",  # set mode, must be either SNP, INDEL or BOTH
#        resources={
#            "hapmap": {"known": False, "training": True, "truth": True, "prior": 15.0},
#            "omni": {"known": False, "training": True, "truth": True, "prior": 12.0},
#            "g1k": {"known": False, "training": True, "truth": True, "prior": 10.0},
#            "dbsnp": {"known": False, "training": False, "truth": False, "prior": 2.0},
#        },
#        annotation=["MQ", "QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR", "DP"],
#        extra=""
#    threads: config['threads']
#    resources:
#        mem_mb=1024,
#    wrapper:
#        "v1.31.1/bio/gatk/variantrecalibrator"




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





