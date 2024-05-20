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
      "Calling long indels using transindel on sample:{wildcards.sample} with group:{wildcards.group}"
    log:
        "logs/{sample}/transindel/{seqtype}_{group}_call.log"
    conda:
        "../envs/transindel.yml"
    params:
        mapq=config['mapq']
    shell:
        """
          python workflow/scripts/transindel/transIndel_call.py \
          -i {input.bam} \
          -l 10 \
          -o results/{wildcards.sample}/{wildcards.seqtype}/indel/transindel/{wildcards.group}_call \
          -m {params} >> {log} 2>&1
        """

# resove alleles and remove PCR slippage
rule long_indel_augment:
    input:
        "results/{sample}/{seqtype}/indel/transindel/{group}_call.indel.vcf"
    output:
        "results/{sample}/{seqtype}/indel/transindel/{group}_long.indels_augmented.vcf"
    message:
      "Augment long indels with group and source information and resolving alleles and removing PCR slippage using transindel on sample:{wildcards.group} with replicate:{wildcards.group}"
    log:
        "logs/{sample}/transindel/{seqtype}_{group}_sliprem.log"
    conda:
        "../envs/manipulate_vcf.yml"
    shell:
        """
            python3 workflow/scripts/add_infos_to_vcf.py \
            {input} \
            long_indel \
            {wildcards.group} \
            {output}_infos >> {log} 2>&1
            
            python3 workflow/scripts/add_contigs_to_vcf.py \
            {output}_infos \
            {output}_contigs \
            resources/refs/genome.fasta >> {log} 2>&1
            
            python3 workflow/scripts/slippage_removal.py \
            resources/refs/genome.fasta \
            {output}_contigs \
            {output} >> {log} 2>&1

            rm -r {output}_contigs {output}_infos

        """

rule longindel_sort_and_compress:
  input:
    "results/{sample}/{seqtype}/indel/transindel/{group}_long.indels_augmented.vcf"
  output:
    "results/{sample}/{seqtype}/indel/transindel/{group}_long.indels.vcf.gz"
  message:
    "Sorting and compressing long indels on sample:{wildcards.sample}"
  log:
    "logs/{sample}/transindel/{seqtype}_{group}_longindel_sort.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools sort {input} -o - | bcftools view -O z -o {output} >> {log} 2>&1
    """

rule combine_longindels:
  input:
    get_longindels
  output:
    "results/{sample}/variants/long.indels.vcf.gz"
  message:
    "Combining long indels on sample:{wildcards.sample}"
  log:
    "logs/{sample}/transindel/combine_longindels.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools concat --naive -O z {input} -o - | bcftools sort -O z -o {output} >> {log} 2>&1
    """

####### MUTECT2 ######

checkpoint split_bam_detect_short_indels_m2:
  input:
    bam="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.baserecal.bam",
    idx="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.baserecal.bam.bai"
  output:
    directory("results/{sample}/{seqtype}/indel/mutect2/{group}_baserecal_split")
  message:
    "Splitting bam file for somatic SNV/Indel detection with Mutect2 on recalibrated data on sample:{wildcards.sample} with group:{wildcards.group}"
  log:
    "logs/{sample}/indel/gatk/mutect2/split_recal_{seqtype}_{group}.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      python workflow/scripts/split_bam_by_chr.py \
          {input.bam} {output} >> {log} 2>&1
    """

rule index_split_bam_detect_short_indels_m2:
  input:
    bam="results/{sample}/{seqtype}/indel/mutect2/{group}_baserecal_split/{chr}.bam"
  output:
    idx="results/{sample}/{seqtype}/indel/mutect2/{group}_baserecal_split/{chr}.bam.bai"
  message:
    "Indexing splitted bam files for somatic SNV/Indel detection with Mutect2 on recalibrated data on sample:{wildcards.sample} with group:{wildcards.group}"
  log:
    "logs/{sample}/indel/gatk/mutect2/index_recal_{seqtype}_{group}_{chr}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      samtools index {input.bam} >> {log} 2>&1
    """

rule detect_short_indels_m2:
  input:
    map="results/{sample}/{seqtype}/indel/mutect2/{group}_baserecal_split/{chr}.bam",
    idx="results/{sample}/{seqtype}/indel/mutect2/{group}_baserecal_split/{chr}.bam.bai",
    fasta="resources/refs/genome.fasta"
  output:
    vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}.vcf"
  message:
    "Detection of somatic SNVs/Indels with Mutect2 on sample:{wildcards.sample} with group:{wildcards.group} on chromosome {wildcards.chr}"
  log:
    "logs/{sample}/indel/gatk/mutect2/m2_{seqtype}_{group}_{chr}.log",
  threads: 4
  resources:
    mem_mb=1024
  wrapper:
    "v1.31.1/bio/gatk/mutect"

rule filter_short_indels_m2:
  input:
    vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}.vcf",
    bam="results/{sample}/{seqtype}/indel/mutect2/{group}_baserecal_split/{chr}.bam",
    idx="results/{sample}/{seqtype}/indel/mutect2/{group}_baserecal_split/{chr}.bam.bai",
    ref="resources/refs/genome.fasta"
  output:
    vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf"
  message:
    "Filtering somatic SNVs/Indels with FilterMutectCalls on sample:{wildcards.sample} with group:{wildcards.group} on chromosome {wildcards.chr}"
  log:
    "logs/{sample}/indel/gatk/filtermutect/filter_{seqtype}_{group}_{chr}.log"
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

rule sort_short_indels_m2:
  input:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf"
  output:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf.gz",
  message:
    "Sorting vcf file from somatic variant calling (mutect2) on recalibrated data on sample:{wildcards.sample} with group:{wildcards.group} on chromosome {wildcards.chr}"
  log:
    "logs/{sample}/gatk/mutect2/sort_{seqtype}_{group}_{chr}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools sort {input} -o - | bcftools view -O z -o {output} >> {log} 2>&1
    """

rule index_short_indels_m2:
  input:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf.gz",
  output:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf.gz.tbi"
  message:
    "Indexing vcf file from somatic variant valling (mutect2) first round on recalibrated data on sample:{wildcards.sample} with group:{wildcards.group} on chromosome {wildcards.chr}"
  log:
     "logs/{sample}/gatk/mutect2/index_{seqtype}_{group}_{chr}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools index -t {input} >> {log} 2>&1
    """

rule merge_short_indels_m2:
  input:
    vcf=aggregate_vcf_mutect2,
    idx=aggregate_idx_mutect2
  output:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz"
  message:
    "Merging vcf files from first round of variant calling (htcaller) on original, unrecalibrated data on sample:{wildcards.sample} with group:{wildcards.group}"
  log:
    "logs/{sample}/gatk/mutect2/merge_{seqtype}_{group}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools concat -O z -a {input.vcf} -o {output} >> {log} 2>&1
    """

######### POST-PROCESSING ########

rule index_merged_short_indels_m2:
  input:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz"
  output:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz.tbi"
  message:
    "Indexing vcf file from somatic variant valling (mutect2) on merged recalibrated data on sample:{wildcards.sample} with group:{wildcards.group}"
  log:
     "logs/{sample}/gatk/mutect2/index_merged_{seqtype}_{group}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools index -t {input} >> {log} 2>&1
    """
 
rule select_short_indels_m2:
    input:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz",
        idx="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz.tbi",
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

rule augment_short_indels_m2:
  input:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf"
  output:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels_augmented.vcf"
  message:
    "Combining somatic short indels detected by Mutect2 on sample:{wildcards.sample}"
  log: 
    "logs/{sample}/combine/{seqtype}/{group}_combine_short.indels.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/add_infos_to_vcf.py \
          {input} \
          short_indel \
          {wildcards.group} \
          {output} >> {log} 2>&1
    """

rule sort_aug_short_indels_m2:
  input:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels_augmented.vcf"
  output:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf.gz"
  message:
    "Sorting and compressing short indels on sample:{wildcards.sample}"
  log:
    "logs/{sample}/transindel/{seqtype}_{group}_shortindel_sort.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools sort {input} -o - | bcftools view -O z -o {output} >> {log} 2>&1
    """

rule combine_aug_short_indels_m2:
  input:
    get_shortindels
  output:
    "results/{sample}/variants/somatic.short.indels.vcf.gz"
  message:
    "Combining long indels on sample:{wildcards.sample}"
  log:
    "logs/{sample}/transindel/combine_longindels.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools concat --naive -O z {input} -o - | bcftools sort -O z -o {output} >> {log} 2>&1
    """
            
rule select_SNVs_m2:
  input:
    vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz",
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

rule augment_somatic_SNVs_m2:
  input:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf"
  output:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs_augmented.vcf"
  message:
    "Combining somatic SNVs detected by Mutect2 on sample:{wildcards.sample}"
  log: 
    "logs/{sample}/indel/{seqtype}/{group}_combine_somatic_SNVs.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/add_infos_to_vcf.py \
          {input} \
          snv \
          {wildcards.group} \
          {output} >> {log} 2>&1
    """

rule sort_somatic_SNVs_m2:
  input:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs_augmented.vcf"
  output:
    "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf.gz"
  message:
    "Sorting and compressing short indels on sample:{wildcards.sample}"
  log:
    "logs/{sample}/transindel/{seqtype}_{group}_SNVs_sort.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools sort {input} -o - | bcftools view -O z -o {output} >> {log} 2>&1
    """

rule combine_somatic_SNVs_m2:
  input:
    get_snvs
  output:
    "results/{sample}/variants/somatic.snvs.vcf.gz"
  message:
    "Combining long indels on sample:{wildcards.sample}"
  log:
    "logs/{sample}/transindel/combine_somatic_SNVs.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools concat --naive -O z {input} -o - | bcftools sort -O z -o {output} >> {log} 2>&1
    """

