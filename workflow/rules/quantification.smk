rule countfeatures_dnaseq:
  input:
    # list of sam or bam files
    sample = "results/{sample}/dnaseq/align/{group}_final_BWA.bam",
    annotation_file = "resources/refs/genome.gtf"
  output:
    table = "results/{sample}/dnaseq/quantification/{group}_counts.txt",
    summary = "results/{sample}/dnaseq/quantification/{group}_counts.txt.summary",
  log:
    "logs/{sample}/featurecounts/dnaseq_{group}.log",
  threads: 2
  conda:
    "../envs/subread.yml"
  shell:
    """
      featureCounts \
          -F GTF \
          -a {input.annotation_file} \
          -t gene \
          -g gene_id \
          --fracOverlap 0.2 \
          -Q {config[mapq]} \
          -T {threads} \
          -o {output.table} {input.sample} > {log} 2>&1
    """

# TODO: add support for PE reads (BWA stores reads as single-end?)

rule countfeatures_rnaseq:
  input:
    # list of sam or bam files
    sample = "results/{sample}/rnaseq/align/{group}_final_STAR.bam",
    annotation_file = "resources/refs/genome.gtf"
  output:
    table = "results/{sample}/rnaseq/quantification/{group}_counts.txt",
    summary = "results/{sample}/rnaseq/quantification/{group}_counts.txt.summary",
  log:
    "logs/{sample}/featurecounts/rnaseq_{group}.log",
  threads: 2
  conda:
    "../envs/subread.yml"
  shell:
    """
    if [ "{config[data][rnaseq_readtype]}" == "SE" ]; then
      featureCounts \
          -F GTF \
          -a {input.annotation_file} \
          -t gene \
          -g gene_id \
          --fracOverlap 0.2 \
          -Q {config[mapq]} \
          -T {threads} \
          -o {output.table} {input.sample} > {log} 2>&1
    elif [ "{config[data][rnaseq_readtype]}" == "PE" ]; then
      featureCounts \
          -p \
          -F GTF \
          -a {input.annotation_file} \
          -t gene \
          -g gene_id \
          --fracOverlap 0.2 \
          -Q {config[mapq]} \
          -T {threads} \
          -o {output.table} {input.sample} > {log} 2>&1
    fi
    """

# merges the count tables from all samples into single table (calculates TPM)
rule merge_countfeatures:
  input:
    get_counts
  output:
    "results/{sample}/quantification/allcounts.txt"
  log:
    "logs/{sample}/quantification/merge_counttables.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      python workflow/scripts/quantification/merge_counttables.py \
          -i '{input}' \
          -n TPM \
          -o {output} > {log} 2>&1
    """

