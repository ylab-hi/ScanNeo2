### align reads to genome using STAR (when reads are in FASTQ)
if config['data']['rnaseq_filetype'] == '.fastq' or config['data']['rnaseq_filetype'] == '.fq':
  rule star_fq_paired_end:
    input:
      unpack(get_star_input),
      idx = "resources/refs/star/",
    output:
      aln="results/{sample}/rnaseq/align/{group}_aligned_STAR.bam",
      log="results/{sample}/rnaseq/align/{group}_aligned_STAR.log",
      sj="results/{sample}/rnaseq/align/{group}_aligned_STAR.tab",
    message:
      "Aligning reads from {wildcards.group} to genome using STAR"
    log:
      "logs/{sample}/star_align/rnaseq_{group}.log"
    params:
      extra=lambda wildcards: f"""--outSAMtype BAM Unsorted \
          --genomeSAindexNbases 10 \
          --outSAMattributes RG \
          --outSAMattrRGline ID:{wildcards.group} \
          --outFilterMultimapNmax 50 \
          --peOverlapNbasesMin 20 \
          --alignSplicedMateMapLminOverLmate 0.5 \
          --alignSJstitchMismatchNmax 5 -1 5 5 \
          --chimOutType WithinBAM HardClip \
          --chimSegmentMin 20 \
          --chimJunctionOverhangMin 10 \
          --chimScoreDropMax 30 \
          --chimScoreJunctionNonGTAG 0 \
          --chimScoreSeparation 1 \
          --chimSegmentReadGapMax 3 \
          --chimMultimapNmax 50"""
    threads: config['threads']
    wrapper:
        "v2.2.1/bio/star/align"

### align reads to genome using STAR (when reads are in BAM - no preprocessing performed)
if config['data']['rnaseq_filetype'] == '.bam':
  checkpoint split_bamfile_RG:
    input:
      unpack(get_star_input),
    output:
      directory("results/{sample}/rnaseq/reads/{group}/bam/")
    conda:
      "../envs/samtools.yml"
    log:
      "logs/{sample}/samtools/split/{group}.logs"
    threads: 10
    shell:
      """
        mkdir -p {output}
        samtools split -@ {threads} \
        -u {output}/noRG.bam \
        -h {input} -f {output}/%!.%. {input}
      """

  rule bamfile_RG_to_fastq:
    input:
      "results/{sample}/rnaseq/reads/{group}/bam/{rg}.bam"
    output:
      "results/{sample}/rnaseq/reads/{group}/fastq/{rg}.fastq"
    message:
      "Converting group:{wildcards.group} BAM file to FASTQ for readgroup:{wildcards.rg}"
    conda:
      "../envs/samtools.yml"
    log:
      "logs/{sample}/samtools/bam2fastq/{group}_{rg}.log"
    threads: config['threads']
    shell:
      """
        samtools collate -Oun128 -@ {threads} {input} \
            | samtools fastq -OT RG -@ {threads} - | gzip -c - > {output}
      """

  rule star_align_bamfile:
    input:
      fq1 ="results/{sample}/rnaseq/reads/{group}/fastq/{rg}.fastq",
      idx ="resources/refs/star/",
    output:
      aln="results/{sample}/rnaseq/align/{group}/{rg}.bam",
      log="results/{sample}/rnaseq/align/{group}/{rg}.log",
      sj="results/{sample}/rnaseq/align/{group}/{rg}.tab"
    log:
      "logs/{sample}/star_align/bam/{group}_{rg}.log"
    params:
      extra=lambda wildcards: f"""--outSAMtype BAM Unsorted --genomeSAindexNbases 10 \
        --readFilesCommand zcat \
        --outSAMattributes RG --outSAMattrRGline ID:{wildcards.rg} \
        --outFilterMultimapNmax 50 --peOverlapNbasesMin 20 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimOutType WithinBAM HardClip --chimSegmentMin 20 \
        --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 --chimMultimapNmax 50"""
    threads: config['threads']
    wrapper:
      "v1.26.0/bio/star/align"
      
  rule merge_alignment_results:
    input:
      aggregate_aligned_rg
    output:
      "results/{sample}/rnaseq/align/{group}_aligned_STAR.bam",
    log:
      "logs/{sample}/samtools/merge/{group}.log",
    params:
      extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
      "v1.32.1/bio/samtools/merge"


# post-processingn and realignment
rule postproc:
  input:
    "results/{sample}/rnaseq/align/{group}_aligned_STAR.bam"
  output:
    "results/{sample}/rnaseq/align/{group}_final_STAR.bam"
  conda:
    "../envs/samtools.yml"
  log:
    "logs/{sample}/postproc/rnaseq_{group}.log"
  threads: 6
  shell:
    """
      samtools view -bh -F 4 --min-MQ {config[mapq]} {input} -o - \
      | samtools sort -n -@ {threads} -m1g -O bam - -o - \
      | samtools fixmate -pcmu -O bam -@ {threads} - - \
      | samtools sort -@ {threads} -m1g -O bam - -o - \
      | samtools markdup -r -@ {threads} - {output} > {log} 2>&1 
      samtools index {output}
    """

rule postproc_bam_index:
  input:
    "results/{sample}/rnaseq/align/{group}_final_STAR.bam"
  output:
    "results/{sample}/rnaseq/align/{group}_final_STAR.bam.bai"
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
    bam="results/{sample}/rnaseq/align/{group}_final_STAR.bam",
    rg="results/{sample}/rnaseq/reads/{group}_readgroups.txt"
  output:
    "results/{sample}/rnaseq/align/{group}_final_BWA.bam"
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

### workflow when aligning paired-end fastq files for DNAseq
rule bwa_align_dnaseq:
  input:
    get_dna_align_input
  output:
    "results/{sample}/dnaseq/align/{group}_aligned_BWA.bam"
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
      -r LB:{wildcards.sample} -r PL:ILLUMINA -r PU:{wildcards.group} - - \
      | samtools sort -@ 6 -n -m1g - -o {output} > {log} 2>&1
    """

rule dnaseq_postproc:
  input:
    "results/{sample}/dnaseq/align/{group}_aligned_BWA.bam"
  output:
    "results/{sample}/dnaseq/align/{group}_final_BWA.bam"
  log:
    "logs/{sample}/postproc/dnaseq_{group}.log"
  conda:
    "../envs/samtools.yml"
  params:
    extra=""
  threads: config['threads']
  shell:
    """
      samtools fixmate -pcmu -O bam -@ 6 {input} - \
          | samtools sort -m1g -O bam - -o - \
          | samtools markdup -r -@ 6 - {output} > {log} 2>&1
      samtools index {output}
    """


