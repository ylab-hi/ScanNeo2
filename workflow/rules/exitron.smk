rule download_genome:
  output:
    "resources/refs/hg38.fa",
    "resources/refs/gencode.v37.annotation.gtf"
  log:
    "logs/prepare_cds.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      curl -L https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz \
          | gzip -d > resources/refs/hg38.fa 
      curl -L  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz \
          | gzip -d > resources/refs/gencode.v37.annotation.gtf
    """
      

rule prepare_cds:
  output:
    "resources/refs/CDS.bed"
  log:
    "logs/prepare_cds.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      cat resources/refs/gencode.v37.annotation.gtf \
          | awk 'OFS="\\t" {{if ($3=="CDS") {{print $1,$4-1,$5,$10,$16,$7}}}}' \
          | tr -d '";' > {output}
    """

rule prepare_scanexitron_config:
  output:
    "resources/scanexitron_config.ini"
  log:
    "logs/prepare_scanexitron_config.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      python3 workflow/scripts/prep_scanexitron_config.py \
          resources/refs/hg38.fa \
          resources/refs/gencode.v37.annotation.gtf \
          resources/refs/CDS.bed \
          {output} > {log}
    """

rule scanexitron:
    input: 
        bam = "results/{sample}/rnaseq/align/{group}_final_STAR.bam",
        fasta = "resources/refs/hg38.fa",
        gtf = "resources/refs/gencode.v37.annotation.gtf",
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
    "results/{sample}/rnaseq/exitron/{group}_exitron.vcf"
  log:
    "logs/exitron2vcf_{sample}_{group}.log"
  conda:
    "../envs/scanexitron.yml"
  shell:
    """
      python workflow/scripts/exitron2vcf2.py \
        {input} \
        {output} \
        resources/refs/hg38.fa > {log}
    """

rule combine_exitrons:
  input:
    get_exitrons,
  output:
    "results/{sample}/rnaseq/exitron/exitrons_combined.vcf"
  message:
    "Combining exitrons on sample:{wildcards.sample}"
  log:
    "logs/{sample}/scanexitron/combine_groups.log"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/combine_vcf.py '{input}' exitron {output} > {log} 2>&1
    """

rule sort_exitrons:
  input:
    "results/{sample}/rnaseq/exitron/exitrons_combined.vcf",
  output:
    "results/{sample}/variants/exitrons.vcf"
  log:
    "logs/{sample}/scanexitron/sort_exitrons.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bcftools sort {input} -o {output}
    """

