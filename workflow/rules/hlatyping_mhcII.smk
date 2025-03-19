######### MHC-II HLA GENOTYPING ###########
rule filter_reads_mhcII_SE:
  input:
    sample=["results/{sample}/hla/reads/{group}_{nartype}_SE.fq"], 
    idx=multiext(
      "resources/hla/bowtie2_index",
      ".1.bt2",
      ".2.bt2",
      ".3.bt2",
      ".4.bt2",
      ".rev.1.bt2",
      ".rev.2.bt2",
    ),
  output:
    "results/{sample}/hla/mhc-II/reads/{group}_{nartype}_flt_SE.bam"
  log:
    "logs/{sample}/genotyping/reads_filtering_mhc-II_{group}_{nartype}.log"
  params:
    extra="",  # optional parameters
  threads: config['threads']  # Use at least two threads
  wrapper:
    "v2.11.1/bio/bowtie2/align"

rule bam2fastq_reads_mhcII_SE:
  input:
    "results/{sample}/hla/mhc-II/reads/{group}_{nartype}_flt_SE.bam"
  output:
    "results/{sample}/hla/mhc-II/reads/{group}_{nartype}_flt_SE.fq"
  log:
    "logs/{sample}/genotyping/{group}_{nartype}_hla_sam.log"
  conda:
    "../envs/samtools.yml"
  threads: 1
  shell:
    """
      samtools fastq -F 4 {input} > {output} 2> {log}
    """

rule filter_reads_mhcII_PE:
  input:
    sample=get_input_filter_reads_mhcII_PE,
    idx=multiext(
      "resources/hla/bowtie2_index",
      ".1.bt2",
      ".2.bt2",
      ".3.bt2",
      ".4.bt2",
      ".rev.1.bt2",
      ".rev.2.bt2",
    ),
  output:
    "results/{sample}/hla/mhc-II/reads/{group}_{nartype}_flt_PE.bam"
  log:
    "logs/{sample}/genotyping/reads_filtering_mhc-II_{group}_{nartype}.log"
  params:
    extra="",  # optional parameters
  threads: config['threads']  # Use at least two threads
  wrapper:
    "v2.11.1/bio/bowtie2/align"

# this rules create the input files for HLA-HD (needs to be PE)
rule finalize_reads_mhcII:
  input:
    get_input_hlatyping_mhcII
  output:
    fwd="results/{sample}/hla/mhc-II/reads/{group}_{nartype}_final_R1.fq",
    rev="results/{sample}/hla/mhc-II/reads/{group}_{nartype}_final_R2.fq"
  log:
    "logs/{sample}/genotyping/finalize_reads_mhcII_{group}_{nartype}.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python workflow/scripts/finalize_mhcII_input.py \
          {input} \
          {output.fwd} \
          {output.rev}
    """


rule hlatyping_mhcII:
  input:
    fwd="results/{sample}/hla/mhc-II/reads/{group}_{nartype}_final_R1.fq",
    rev="results/{sample}/hla/mhc-II/reads/{group}_{nartype}_final_R2.fq"
  output:
    "results/{sample}/hla/mhc-II/genotyping/{group}_{nartype}/result/{group}_{nartype}_final.result.txt"
  log:
    "logs/{sample}/hla/{group}_{nartype}_hlahd.log"
  conda:
    "../envs/hlahd.yml"
  params:
    freqdata=f"""-f {config['hlatyping']['freqdata']}""",
    split=f"""{config['hlatyping']['split']}""",
    dic=f"""{config['hlatyping']['dict']}"""
  threads: config['threads']
  shell:
    """
      hlahd.sh \
          -t {threads} \
          -m 100 \
          -c 0.95 \
          {params.freqdata} \
          {input.fwd} {input.rev} \
          {params.split} {params.dic} \
          {wildcards.group}_{wildcards.nartype} \
          results/{wildcards.sample}/hla/mhc-II/genotyping/
    """

rule merge_predicted_mhcII_allels:
  input:
    get_predicted_mhcII_alleles
  output:
    "results/{sample}/hla/mhc-II/genotyping/mhc-II.tsv",
  message:
    "Merging HLA alleles from different sources"
  log:
    "logs/{sample}/genotyping/merge_predicted_mhc-II.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python workflow/scripts/genotyping/merge_predicted_mhcII.py \
          '{input}' {output} > {log} 2>&1
    """


# combine predicted (and user-defined) alleles
rule combine_all_mhcII_alleles:
  input:
    get_all_mhcII_alleles
  output:
    "results/{sample}/hla/mhc-II.tsv"
  message:
    "Combining HLA mhc-II alleles from different sources"
  log:
    "logs/{sample}/genotyping/combine_all_mhc-II.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python workflow/scripts/genotyping/combine_all_alleles.py \
          '{input}' mhc-II {output} > {log} 2>&1
    """
