rule prepare_cds:
  input:
    "resources/refs/genome.gtf"
  output:
    "resources/refs/CDS.bed"
  log:
    "logs/prepare_cds.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      cat resources/refs/genome.gtf \
          | awk 'OFS="\\t" {{if ($3=="CDS") {{print $1,$4-1,$5,$10,$16,$7}}}}' \
          | tr -d '";' > {output}
    """

rule prepare_scanexitron_config:
  input:
    genome="resources/refs/genome.fasta",
    annotation="resources/refs/genome.gtf",
    cds="resources/refs/CDS.bed"
  output:
    "resources/scanexitron_config.ini"
  log:
    "logs/prepare_scanexitron_config.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      python3 workflow/scripts/prep_scanexitron_config.py \
          {input.genome} \
          {input.annotation} \
          {input.cds} \
          {output} > {log}
    """
      
rule scanexitron:
    input: 
        bam = "results/{sample}/rnaseq/align/{group}_final_STAR.bam",
        idx = "results/{sample}/rnaseq/align/{group}_final_STAR.bam.bai",
        fasta = "resources/refs/genome.fasta",
        gtf = "resources/refs/genome.gtf",
        cds = "resources/refs/CDS.bed",
        config="resources/scanexitron_config.ini"
    output:
        "results/{sample}/rnaseq/exitron/{group}.exitron",
    message:
      "Detect exitrons on sample:{wildcards.sample} of group:{wildcards.group}"
    log:
        "logs/scanexitron_{sample}_{group}.log"
    threads:
      config["threads"]
    container:
      "docker://yanglabinfo/scanneo2-scanexitron"
    params:
      mapq = config['mapq'],
      ao = config['exitronsplicing']['ao'],
      pso = config['exitronsplicing']['pso']
    shell:
      """
        python3 workflow/scripts/scanexitron/ScanExitron.py \
          -t {threads} \
          --mapq {params.mapq} \
          --ao {params.ao} \
          --pso {params.pso} \
          -c ../../../{input.config} \
          -i {input.bam} \
          -r hg38
        mv {wildcards.group}_final_STAR.exitron {output}
        mv {wildcards.group}_final_STAR* results/{wildcards.sample}/rnaseq/exitron/
      """

rule exitron_to_vcf:
  input:
    "results/{sample}/rnaseq/exitron/{group}.exitron"
  output:
    "results/{sample}/rnaseq/exitron/{group}_exitrons.vcf"
  log:
    "logs/exitron2vcf_{sample}_{group}.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/exitron2vcf.py \
        {input} {output} \
        resources/refs/genome.fasta > {log}
    """

rule exitron_augment:
  input:
    "results/{sample}/rnaseq/exitron/{group}_exitrons.vcf"
  output:
    "results/{sample}/rnaseq/exitron/{group}_exitrons_augmented.vcf"
  message:
    "Augmenting exitrons on sample:{wildcards.sample} of group:{wildcards.group}"
  log:
    "logs/exitron_augment_{sample}_{group}.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/add_infos_to_vcf.py {input} exitron {output} > {log} 2>&1
    """

rule sort_exitron:
  input:
    "results/{sample}/rnaseq/exitron/{group}_exitrons_augmented.vcf"
  output:
    "results/{sample}/rnaseq/exitron/{group}_exitrons.vcf.gz"
  message:
    "Sorting and compressing exitrons on sample:{wildcards.sample} of group:{wildcards.group}"
  log:
    "logs/exitron_sort_{sample}_{group}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools sort {input} -o - | bcftools view -O z -o {output} > {log} 2>&1
    """

rule combine_exitrons:
  input:
    get_exitrons
  output:
    "results/{sample}/variants/exitrons.vcf.gz"
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

