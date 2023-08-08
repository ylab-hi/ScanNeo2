rule postproc:
  input:
    "results/{sample}/rnaseq/align/{group}_aligned.bam"
  output:
    "results/{sample}/rnaseq/align/{group}_ready.bam"
  conda:
    "../envs/samtools.yml"
  log:
    "logs/{sample}/postproc/rnaseq_{group}.log"
  threads: 6
  shell:
    """
      samtools index {input}
      samtools view -bh -F 4 --min-MQ {config[mapq]} {input} -o - \
      | samtools sort -n -@ {threads} -m1g -O bam - -o - \
      | samtools fixmate -pcmu -O bam -@ {threads} - - \
      | samtools sort -@ {threads} -m1g -O bam - -o - \
      | samtools markdup -r -@ {threads} - {output} > {log} 2>&1 
      samtools index {output}
    """

rule postproc_bam_index:
  input:
    "results/{sample}/rnaseq/align/{group}_ready.bam"
  output:
    "results/{sample}/rnaseq/align/{group}_ready.bam.bai"
  conda:
    "../envs/samtools.yml"
  log:
    "logs/{sample}/postproc/index/rnaseq_{group}.log"
  shell:
    """
      samtools index {input} > {log} 2>&1
    """


## retrieve readgroups from bam file
rule get_readgroups:
  input:
    get_readgroups_input
  output:
        "results/{sample}/{seqtype}/reads/{group}_readgroups.txt"
  conda:
      "../envs/basic.yml"
  log:
        "logs/{sample}/get_readgroups/{seqtype}_{group}.log"
  shell:
    """
          python workflow/scripts/get_readgroups.py '{input}' \
          {output} > {log} 2>&1
      """

rule realign:
  input:
    bam="results/{sample}/rnaseq/align/{group}_ready.bam",
    rg="results/{sample}/rnaseq/reads/{group}_readgroups.txt"
  output:
    "results/{sample}/rnaseq/align/{group}_realigned.bam"
  conda:
    "../envs/basic.yml"
  log:
    "logs/{sample}/realign/rnaseq_{group}.log"
  threads: config['threads']
  shell:
    """
      samtools collate -Oun128 {input.bam} \
        | samtools fastq -OT RG -@ {threads} - \
        | bwa mem -pt{threads} -CH <(cat {input.rg}) resources/refs/bwa/genome - \
        | samtools sort -@6 -m1g - -o {output} > {log} 2>&1
        samtools index {output}
    """
