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
      #\
          #| gzip -d > resources/refs/gencode.v37.annotation.gtf
      

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
      """

rule exitron_to_vcf:
  input:
    "results/exitron/{sample}/{group}.exitron"
  output:
    "results/exitron/{sample}/{group}.vcf"
  log:
    "logs/exitron2vcf_{sample}_{group}.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      workflow/scripts/exitron2vcf.py \
        resources/refs/genome.fasta \
        '{input}' \
        {output} 2> {log}
    """

#rule combine_exitrons:
    #input:
        #expand("results/exitron/{sample}/res.vcf", sample=config["rnaseq"])
    #output:
        #"results/exitron/all.vcf"
    #shell:
        #"""
        #"""
