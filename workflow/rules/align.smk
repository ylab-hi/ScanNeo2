### align reads to genome using STAR (when reads are in FASTQ)
if config['data']['rnaseq_filetype'] == '.fastq' or config['data']['rnaseq_filetype'] == '.fq':
  rule star_align_fastq:
    input:
      unpack(get_star_input),
      faidx = "resources/refs/genome.fasta.fai",
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
          --outSAMattributes RG HI \
          --outSAMattrRGline ID:{wildcards.group} \
          --outFilterMultimapNmax 50 \
          --peOverlapNbasesMin 15 \
          --alignSplicedMateMapLminOverLmate 0.5 \
          --alignSJstitchMismatchNmax 5 -1 5 5 \
          --chimOutType WithinBAM HardClip \
          --chimSegmentMin {config["align"]["chimSegmentMin"]} \
          --chimJunctionOverhangMin {config["align"]["chimJunctionOverhangMin"]} \
          --chimScoreDropMax {config["align"]["chimScoreDropMax"]} \
          --chimScoreMin {config["align"]["chimScoreMin"]} \
          --chimScoreJunctionNonGTAG 0 \
          --chimScoreSeparation {config["align"]["chimScoreSeparation"]} \
          --chimSegmentReadGapMax 3 \
          --chimMultimapNmax 50 \
          --outSAMstrandField intronMotif"""
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
      faidx = "resources/refs/genome.fasta.fai",
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
        --outSAMattributes RG HI --outSAMattrRGline ID:{wildcards.rg} \
        --outFilterMultimapNmax 50 \
        --peOverlapNbasesMin 15 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimOutType WithinBAM HardClip \
        --chimSegmentMin {config["align"]["chimSegmentMin"]} \
        --chimJunctionOverhangMin {config["align"]["chimJunctionOverhangMin"]} \
        --chimScoreDropMax {config["align"]["chimScoreDropMax"]} \
        --chimScoreMin {config["align"]["chimScoreMin"]} \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation {config["align"]["chimScoreSeparation"]} \
        --chimSegmentReadGapMax 3 --chimMultimapNmax 50 \
        --outSAMstrandField intronMotif"""
    threads: config['threads']
    wrapper:
      "v1.26.0/bio/star/align"
      
  rule merge_alignment_results:
    input:
      aggregate_aligned_rg
    output:
      "results/{sample}/rnaseq/align/{group}_aligned_STAR.bam"
    log:
      "logs/{sample}/samtools/merge/{group}.log"
    params:
      extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
      "v1.32.1/bio/samtools/merge"

# post-processingn and realignment
rule rnaseq_postproc_fixmate:
  input:
    "results/{sample}/rnaseq/align/{group}_aligned_STAR.bam"
  output:
    "results/{sample}/rnaseq/align/{group}_fixmate_STAR.bam"
  conda:
    "../envs/samtools.yml"
  log:
    "logs/{sample}/postproc/rnaseq_{group}_fixmate.log"
  threads: 4
  params:
    mapq="--min-MQ config['mapq']"
  shell:
    """
      samtools view -h -F 4 {params.mapq} {input} -o - \
      | samtools sort -n -@4 -m4g -O SAM - -o - \
      | samtools fixmate -pcmu -O bam -@ {threads} - {output} > {log} 2>&1
    """

# sort and markdup needed to be separated (ensure no core dump for whatever reason) 
rule rnaseq_postproc_markdup:
  input:
    bam="results/{sample}/rnaseq/align/{group}_fixmate_STAR.bam",
    tmp="tmp/"
  output:
    "results/{sample}/rnaseq/align/{group}_final_STAR.bam"
  conda:
    "../envs/samtools.yml"
  log:
    "logs/{sample}/postproc/rnaseq_{group}_markdup.log"
  threads: 4
  resources:
    mem_mb_per_cpu=4000
  shell:
    """
      samtools sort -@4 -m4G -O BAM -T tmp/ {input.bam} \
          -o tmp/rnaseq_fixmate_sorted_{wildcards.sample}_{wildcards.group}.bam > {log} 2>&1
      samtools markdup -r -@4 tmp/rnaseq_fixmate_sorted_{wildcards.sample}_{wildcards.group}.bam \
          {output} > {log} 2>&1 
      rm tmp/rnaseq_fixmate_sorted_{wildcards.sample}_{wildcards.group}.bam
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

# realign RNAseq and align DNAseq 
rule realign:
  input:
    bam=get_readgroups_input,
    rg="results/{sample}/{seqtype}/reads/{group}_readgroups.txt",
    idx = multiext("resources/refs/bwa/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
  output:
    bam="results/{sample}/{seqtype}/align/{group}_final_BWA.bam",
  conda:
    "../envs/realign.yml"
  log:
    "logs/{sample}/realign/{seqtype}_{group}.log"
  threads: config['threads']
  shell:
    """
      samtools collate -Oun128 {input.bam} \
        | samtools fastq -OT RG -@ {threads} - \
        | bwa mem -pt{threads} -CH <(cat {input.rg}) resources/refs/bwa/genome - - \
        | samtools sort -@6 -m1g -o {output} > {log} 2>&1
    """


### workflow when aligning paired-end fastq files for DNAseq
if config['data']['dnaseq_filetype'] in ['.fq','.fastq']:
  rule bwa_align_dnaseq:
    input:
      reads=get_dna_align_input,
      idx = multiext("resources/refs/bwa/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
      "results/{sample}/dnaseq/align/{group}_aligned_BWA.bam"
    log:
      "logs/{sample}/bwa_align/dnaseq_{group}.log"
    conda:
      "../envs/realign.yml"
    params:
      extra=""
    threads: config['threads']
    shell:
      """
        bwa mem -t{threads} resources/refs/bwa/genome \
            -R '@RG\\tID:{wildcards.group}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA' \
            {input.reads} | samtools sort -@ 6 -n -m1g - -o {output} > {log} 2>&1
      """

  rule dnaseq_postproc:
    input:
      aln="results/{sample}/dnaseq/align/{group}_aligned_BWA.bam",
      tmp="tmp/"
    output:
      bam="results/{sample}/dnaseq/align/{group}_final_BWA.bam",
    log:
      "logs/{sample}/postproc/dnaseq_{group}.log"
    conda:
      "../envs/samtools.yml"
    params:
      extra=""
    threads: config['threads']
    shell:
      """
        samtools fixmate -pcmu -O bam -@ 6 {input.aln} - \
            | samtools sort -m1g -O bam -T tmp/ - -o - \
            | samtools markdup -r -@ 6 - {output.bam} > {log} 2>&1
      """
    
rule samtools_index_BWA_final:
    input:
        "results/{sample}/{seqtype}/align/{group}_final_BWA.bam",
    output:
        "results/{sample}/{seqtype}/align/{group}_final_BWA.bam.bai",
    log:
        "logs/samtools_index/realign_{seqtype}_{sample}_{group}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v2.3.0/bio/samtools/index"
