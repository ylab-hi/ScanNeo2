rule countfeatures:
  input:
    sample = get_aligned_reads,
    annotation_file = "resources/refs/genome.gtf"
  output:
    "results/{sample}/{seqtype}/quantification/{group}_counts.txt"
  log:
    "logs/{sample}/featurecounts/{seqtype}_{group}.log"
  threads: 2
  conda:
    "../envs/subread.yml"
  params:
    mapq=f"""{config['mapq']}"""
  shell:
    """
      python workflow/scripts/quantification/featurecounts_wrapper.py \
          {input.sample} {output} {input.annotation_file} \
          {params.mapq} {threads} 
    """

# merges the count tables from all samples into single table (calculates TPM)
rule merge_counttables:
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
