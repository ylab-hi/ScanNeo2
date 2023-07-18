### workflow when aligning paired-end fastq files for RNAseq
checkpoint splitfastq:
  input:
    get_splitfastq_input
  output:
    directory("results/{sample}/rnaseq/reads/{group}/")
  log:
    "logs/{sample}/splitfastq/rnaseq_{group}.log"
  conda:
    "../envs/splitfastq.yml"
  threads: 0
  shell:
    """
      python workflow/scripts/splitfastq.py '{input}' {output} 20000000
    """

rule star_fq_paired_end:
  input:
    fq1 = "results/{sample}/rnaseq/reads/{group}/r1/reads_{i}.fq.gz",
    fq2 = "results/{sample}/rnaseq/reads/{group}/r2/reads_{i}.fq.gz",
    idx = "resources/refs/star/",
  output:
    aln = "results/{sample}/rnaseq/align/{group}/reads_{i}.bam",
    log = "results/{sample}/rnaseq/align/{group}/reads_{i}.log",
    sj = "results/{sample}/rnaseq/align/{group}/reads_{i}.tab"
  log:
    "logs/{sample}/star_align/rnaseq_{group}_{i}.log"
  params:
    extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10 --outSAMattributes RG --outSAMattrRGline ID:{group}"
  threads: config['threads']
  wrapper:
      "v1.26.0/bio/star/align"

rule merge_align:
  input:
    aggregate_align,
  output:
    "results/{sample}/rnaseq/align/{group}_aligned.bam"
  log:
    "logs/{sample}/star_align/rnaseq_{group}.log"
  params:
      #extra="",  
  threads: config['threads']
  wrapper:
    "v1.32.1/bio/samtools/merge"


### workflow when aligning paired-end fastq files for DNAseq
rule bwa_align_dnaseq:
  input:
    get_dna_align_input
  output:
    "results/{sample}/dnaseq/align/{group}_realigned.bam"
  log:
    "logs/{sample}/bwa_align/dnaseq_{group}.log"
  conda:
    "../envs/basic.yml"
  params:
    extra=""
  threads: config['threads']
  shell:
    """
      bwa mem -pt{threads} -C resources/refs/bwa/genome {input} \
      | samtools addreplacerg -r ID:{wildcards.group} -r SM:{wildcards.sample} \
      -r LB:{wildcards.sample} -r PL:ILLUMINA -r PU:{wildcards.group} - \
      | samtools fixmate -pcmu -O bam -@ {threads} - - \
      | samtools sort -@ {threads} -m1g -O bam - -o - \
      | samtools markdup -r -@ {threads} - {output} > {log} 2>&1 
    """
