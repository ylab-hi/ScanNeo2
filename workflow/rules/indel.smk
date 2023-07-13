import os
from snakemake.remote import HTTP

rule detect_long_indel_ti_build:
    input:
        bam = "results/{sample}/rnaseq/align/{replicate}_realigned.bam",
        idx = "results/{sample}/rnaseq/align/{replicate}_realigned.bam.bai"
    output:
        bam="results/{sample}/rnaseq/indel/transindel/{replicate}_build.bam",
        idx="results/{sample}/rnaseq/indel/transindel/{replicate}_build.bam.bai"
    message:
      "Building new BAM file with redefined CIGAR string using transindel build on sample:{wildcards.sample} with replicate:{wildcards.replicate}"
    log:
        "logs/{sample}/transindel/{replicate}_build.log"
    conda:
        "../envs/transindel.yml"
    shell:
        """
        python3 workflow/scripts/transIndel/transIndel_build_RNA.py \
        -i {input.bam} \
        -o {output.bam} \
        -r resources/refs/genome.fasta \
        -g resources/refs/genome.gtf > {log} 2>&1
        samtools index {output.bam} -o {output.idx} >> {log} 2>&1
        """

rule detect_long_indel_ti_call:
    input:
        bam = "results/{sample}/rnaseq/indel/transindel/{replicate}_build.bam",
        bai = "results/{sample}/rnaseq/indel/transindel/{replicate}_build.bam.bai"
    output:
        "results/{sample}/rnaseq/indel/transindel/{replicate}_call.indel.vcf"
    message:
      "Calling short indels using transindel on sample:{wildcards.sample} with replicate:{wildcards.replicate}"
    log:
        "logs/{sample}/transindel/{replicate}_call.log"
    conda:
        "../envs/transindel.yml"
    params:
        mapq=config['mapq']
    shell:
        """
        python workflow/scripts/transIndel/transIndel_call.py \
        -i {input.bam} \
        -l 10 \
        -o results/{wildcards.sample}/rnaseq/indel/transindel/{wildcards.replicate}_call \
        -m {params} > {log} 2>&1
        """

# resove alleles and remove PCR slippage
rule long_indel_slippage_removal:
    input:
        "results/{sample}/rnaseq/indel/transindel/{replicate}_call.indel.vcf"
    output:
        "results/{sample}/rnaseq/indel/transindel/{replicate}_sliprem.vcf"
    message:
      "Resolving alleles and removing PCR slippage using transindel on sample:{wildcards.sample} with replicate:{wildcards.replicate}"
    log:
        "logs/{sample}/transindel/{replicate}_sliprem.log"
    conda:
        "../envs/transindel.yml"
    shell:
        """
            python3 workflow/scripts/slippage_removal.py \
            resources/refs/genome.fasta {input} {output} > {log} 2>&1
        """

# combines the replicates into one vcf
rule combine_longindels:
  input:
    expand("results/{sample}/rnaseq/indel/transindel/{replicate}_sliprem.vcf", 
           sample=config['data']['name'],
           replicate=config['data']['rnaseq'].keys())
  output:
    "results/{sample}/rnaseq/indel/long.indel.vcf" 
  message:
    "Combining long indels from replicates on sample:{wildcards.sample}"
  log: 
    "logs/{sample}/transindel/combine_replicates.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/combine_vcf.py '{input}' {output} > {log} 2>&1
    """

# detects short somatic variants (SNVs and indels) using mutect2
rule detect_short_indels_m2:
    input:
        fasta="resources/refs/genome.fasta",
        map="results/{sample}/rnaseq/align/{replicate}_realigned.bam"
    output:
        vcf="results/{sample}/rnaseq/indel/mutect2/{replicate}_variants.vcf",
        bam="results/{sample}/rnaseq/indel/mutect2/{replicate}_variants.bam",
    message:
      "Detection of somatic SNVs/Indels with Mutect2 on sample:{wildcards.sample} with replicate:{wildcards.replicate}"
    threads: config['threads']
    resources:
        mem_mb=10024,
    params:
        extra="",
    log:
        "logs/{sample}/gatk/mutect2/{replicate}.log",
    wrapper:
        "v1.31.1/bio/gatk/mutect"


# filters short somatic variants (SNVs and indels) using FilterMutectCalls
rule filter_short_indels:
    input:
        vcf="results/{sample}/rnaseq/indel/mutect2/{replicate}_variants.vcf",
        bam="results/{sample}/rnaseq/indel/mutect2/{replicate}_variants.bam",
        ref="resources/refs/genome.fasta",
        # intervals="intervals.bed",
        # contamination="", # from gatk CalculateContamination
        # segmentation="", # from gatk CalculateContamination
        # f1r2="", # from gatk LearnReadOrientationBias
    output:
        vcf="results/{sample}/rnaseq/indel/mutect2/{replicate}_variants.flt.vcf"
    message:
      "Filtering somatic SNVs/Indels with FilterMutectCalls on sample:{wildcards.sample} and replicate:{wildcards.replicate}"
    log:
        "logs/{sample}/gatk/filtermutect/{replicate}.log",
    params:
        extra="--max-alt-allele-count 3",
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/filtermutectcalls"

rule combine_short_indels_m2:
  input:
    expand("results/{sample}/rnaseq/indel/mutect2/{replicate}_variants.flt.vcf", 
           sample=config['data']['name'],
           replicate=config['data']['rnaseq'].keys())
  output:
    "results/{sample}/rnaseq/indel/mutect2/variants.vcf" 
  message:
    "Combining somatic SNVs/Indels with Mutect2 on sample:{wildcards.sample}"
  log: 
    "logs/{sample}/transindel/combine_replicates.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/combine_vcf.py '{input}' {output} > {log} 2>&1
    """

rule select_SNVs_m2:
    input:
        vcf="results/{sample}/rnaseq/indel/mutect2/variants.vcf",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/rnaseq/indel/snvs.vcf"
    message:
      "Selecting somatic SNVs with SelectVariants on sample:{wildcards.sample}"
    log:
        "logs/{sample}/gatk/select/somatic_snvs.log",
    params:
        extra="--select-type-to-include SNP",  # optional filter arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/selectvariants"

rule select_short_indels_m2:
    input:
        vcf="results/{sample}/rnaseq/indel/mutect2/variants.vcf",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/rnaseq/indel/short.indel.vcf"
    message:
      "Selecting short somatic indels with SelectVariants on sample:{wildcards.sample}"
    log:
        "logs/gatk/select/{sample}.indel.log",
    params:
        extra="--select-type-to-include INDEL",  # optional filter arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/selectvariants"

