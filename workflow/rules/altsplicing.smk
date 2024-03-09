rule spladder:
    input:
        bam = "results/{sample}/rnaseq/align/{group}_final_STAR.bam",
        bamidx = "results/{sample}/rnaseq/align/{group}_final_STAR.bam.bai"
    output:
        directory("results/{sample}/rnaseq/altsplicing/spladder/{group}")
    conda: 
        "../envs/spladder.yml"
    log:
        "logs/{sample}/spladder/{group}_build.log"
    params:
      confidence = "--confidence {config[altsplicing][confidence]}",
      iteration = "--iterations {config[altsplicing][iterations]}",
      edgelimit = "--ase-edge-limit {config[altsplicing][edgelimit]}",
    shell:
        """
          spladder build -b {input.bam} \
              -a resources/refs/genome.gtf \
              -o {output} --filter-overlap-exons \
              --no-primary-only --quantify-graph \
              {params.confidence} \
              {params.iteration} \
              {params.edgelimit} \
              --qmode all > {log} 2>&1
        """

rule splicing_to_vcf:
  input:
    "results/{sample}/rnaseq/altsplicing/spladder/{group}"
  output:
    "results/{sample}/rnaseq/altsplicing/spladder/{group}_altsplicing.vcf"
  message:
    "Converting splicing events to VCF format"
  log:
    "logs/{sample}/spladder/{group}_to_vcf.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/altsplc2vcf.py \
          -i {input} -r resources/refs/genome.fasta \
          -g {wildcards.group} -o {output} > {log} 2>&1
    """

rule sort_altsplicing:
  input:
    "results/{sample}/rnaseq/altsplicing/spladder/{group}_altsplicing.vcf"
  output:
    "results/{sample}/rnaseq/altsplicing/spladder/{group}_altsplicing.vcf.gz"
  message:
    "Sorting and compressing splicing events on sample:{wildcards.sample} of group:{wildcards.group}"
  log:
    "logs/{sample}/spladder/{group}_sort.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools sort {input} -o - | bcftools view -O z -o {output} > {log} 2>&1
    """

rule combine_altsplicing:
  input:
    get_altsplicing
  output:
    "results/{sample}/variants/altsplicing.vcf.gz"
  message:
    "Combining exitrons on sample:{wildcards.sample}"
  log:
    "logs/{sample}/exitrons/combine_exitrons.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools concat --naive -O z {input} -o  - | bcftools sort -O z -o {output} > {log} 2>&1
    """
