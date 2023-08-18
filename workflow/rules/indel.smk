import os
from snakemake.remote import HTTP
rule detect_long_indel_ti_build_RNA:
    input:
        bam = "results/{sample}/rnaseq/align/{group}_final_BWA.bam",
        idx = "results/{sample}/rnaseq/align/{group}_final_BWA.bam.bai"
    output:
        bam="results/{sample}/rnaseq/indel/transindel/{group}_build.bam",
        idx="results/{sample}/rnaseq/indel/transindel/{group}_build.bam.bai"
    message:
      "Building new BAM file with redefined CIGAR string using transindel build on sample:{wildcards.sample} with group:{wildcards.group}"
    log:
        "logs/{sample}/transindel/{group}_build_RNA.log"
    conda:
        "../envs/transindel.yml"
    shell:
        """
          python3 workflow/scripts/transindel/transIndel_build_RNA.py \
          -i {input.bam} \
          -o {output.bam} \
          -r resources/refs/genome.fasta \
          -g resources/refs/genome.gtf > {log} 2>&1
          samtools index {output.bam} -o {output.idx} >> {log} 2>&1
        """

rule detect_long_indel_ti_build_DNA:
    input:
        bam = "results/{sample}/dnaseq/align/{group}_final_BWA.bam",
        idx = "results/{sample}/dnaseq/align/{group}_final_BWA.bam.bai"
    output:
        bam="results/{sample}/dnaseq/indel/transindel/{group}_build.bam",
        idx="results/{sample}/dnaseq/indel/transindel/{group}_build.bam.bai"
    message:
      "Building new BAM file with redefined CIGAR string using transindel build on sample:{wildcards.sample} with group:{wildcards.group}"
    log:
        "logs/{sample}/transindel/{group}_build_DNA.log"
    conda:
        "../envs/transindel.yml"
    shell:
        """
          python workflow/scripts/transindel/transIndel_build_DNA.py \
          -i {input.bam} -o {output.bam}  > {log} 2>&1
          samtools index {output.bam} -o {output.idx} >> {log} 2>&1
        """

rule detect_long_indel_ti_call:
    input:
        bam = "results/{sample}/{seqtype}/indel/transindel/{group}_build.bam",
        bai = "results/{sample}/{seqtype}/indel/transindel/{group}_build.bam.bai"
    output:
        "results/{sample}/{seqtype}/indel/transindel/{group}_call.indel.vcf"
    message:
      "Calling short indels using transindel on sample:{wildcards.sample} with group:{wildcards.group}"
    log:
        "logs/{sample}/transindel/{seqtype}_{group}_call.log"
    conda:
        "../envs/transindel.yml"
    params:
        mapq=config['mapq']
    shell:
        """
          python workflow/scripts/transIndel/transIndel_call.py \
          -i {input.bam} \
          -l 10 \
          -o results/{wildcards.sample}/{wildcards.seqtype}/indel/transindel/{wildcards.group}_call \
          -m {params} > {log} 2>&1
        """

# resove alleles and remove PCR slippage
rule long_indel_slippage_removal:
    input:
        "results/{sample}/{seqtype}/indel/transindel/{group}_call.indel.vcf"
    output:
        "results/{sample}/{seqtype}/indel/transindel/{group}_sliprem.vcf"
    message:
      "Resolving alleles and removing PCR slippage using transindel on sample:{wildcards.group} with replicate:{wildcards.group}"
    log:
        "logs/{sample}/transindel/{seqtype}_{group}_sliprem.log"
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
    get_longindels,
  output:
    "results/{sample}/variants/long.indels.vcf" 
  message:
    "Combining long indels from replicates on sample:{wildcards.sample}"
  log: 
    "logs/{sample}/transindel/combine_groups.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/combine_vcf.py '{input}' long_indel {output} > {log} 2>&1
    """


# detects short somatic variants (SNVs and indels) using mutect2
rule detect_short_indels_m2:
    input:
        fasta="resources/refs/genome.fasta",
        map="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.baserecal.bam",
    output:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf",
        bam="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.bam",
    message:
      "Detection of somatic SNVs/Indels with Mutect2 on sample:{wildcards.sample} with group:{wildcards.group}"
    threads: config['threads']
    resources:
        mem_mb=10024,
    params:
        extra="",
    log:
        "logs/{sample}/gatk/mutect2/{seqtype}_{group}.log",
    wrapper:
        "v1.31.1/bio/gatk/mutect"


# filters short somatic variants (SNVs and indels) using FilterMutectCalls
rule filter_short_indels:
    input:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf",
        bam="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.baserecal.bam",
        ref="resources/refs/genome.fasta",
        # intervals="intervals.bed",
        # contamination="", # from gatk CalculateContamination
        # segmentation="", # from gatk CalculateContamination
        # f1r2="", # from gatk LearnReadOrientationBias
    output:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.flt.vcf"
    message:
      "Filtering somatic SNVs/Indels with FilterMutectCalls on sample:{wildcards.sample} and group:{wildcards.group}"
    log:
        "logs/{sample}/gatk/filtermutect/{seqtype}_{group}.log",
    params:
        extra=f"""--max-alt-allele-count 3 \
          --min-median-base-quality {config['basequal']} \
          --min-median-mapping-quality {config['mapq']} \
          --threshold-strategy {config['indel']['strategy']} \
          --f-score-beta {config['indel']['fscorebeta']} \
          --false-discovery-rate {config['indel']['fdr']} \
          --pcr-slippage-rate {config['indel']['sliprate']} \
          --min-slippage-length {config['indel']['sliplen']}""",
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/filtermutectcalls"


rule select_short_indels_m2:
    input:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.flt.vcf",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf"
    message:
      "Selecting short somatic indels with SelectVariants on sample:{wildcards.sample}"
    log:
        "logs/{sample}/gatk/select/{seqtype}_{group}.indel.log",
    params:
        extra="--select-type-to-include INDEL",  # optional filter arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.31.1/bio/gatk/selectvariants"

rule combine_short_indels_m2:
  input:
    get_shortindels,
  output:
    "results/{sample}/variants/somatic.short.indels.vcf"
  message:
    "Combining somatic short indels detected by Mutect2 on sample:{wildcards.sample}"
  log: 
    "logs/{sample}/transindel/combine_replicates.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/combine_vcf.py '{input}' short_indel {output} > {log} 2>&1
    """

rule select_SNVs_m2:
  input:
    vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.flt.vcf",
    ref="resources/refs/genome.fasta",
  output:
    vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf"
  message:
    "Selecting somatic SNVs with SelectVariants on sample:{wildcards.sample}"
  log:
    "logs/{sample}/gatk/select/{seqtype}_{group}_somatic_snvs.log",
  params:
    extra="--select-type-to-include SNP",  # optional filter arguments, see GATK docs
    java_opts="",  # optional
  resources:
    mem_mb=1024,
  wrapper:
    "v1.31.1/bio/gatk/selectvariants"

rule combine_snvs_m2:
  input:
    get_snvs,
  output:
    "results/{sample}/variants/somatic.snvs.vcf"
  message:
    "Combining somatic SNVs with Mutect2 on sample:{wildcards.sample}"
  log: 
    "logs/{sample}/transindel/combine_replicates.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/combine_vcf.py '{input}' snv {output} > {log} 2>&1
    """

