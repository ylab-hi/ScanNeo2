rule download_vep_plugins:
  output:
    directory("resources/vep/plugins")
  message:
    "Downloading VEP plugins"
  log:
    "logs/vep/plugins.log"
  params:
    release=110
  wrapper:
    "v1.31.1/bio/vep/plugins"

rule download_wildtype_plugin:
  output:
    "resources/vep/plugins/Wildtype.pm"
  message:
    Downloading VEP Wildtype plugin"
  log:
    "logs/vep/plugins_wildtype.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      curl -L -o {output} https://raw.githubusercontent.com/griffithlab/pVAC-Seq/master/pvacseq/VEP_plugins/Wildtype.pm
    """

rule download_vep_cache:
  output:
    directory("resources/vep/cache")
  conda:
    "../envs/basic.yml"
  log:
    "logs/vep/cache.log"
  shell:
    """
      mkdir -p {output}
      curl -L -o - https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz 
      | tar -xz -C resources/vep/cache
    """

rule index_variants:
  input:
    "results/{sample}/variants/{vartype}.vcf"
  output:
    "results/{sample}/variants/{vartype}.vcf.gz",
    "results/{sample}/variants/{vartype}.vcf.gz.tbi"
  log:
    "logs/indexvcf/{sample}_{vartype}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      bgzip {input}
      tabix {input}.gz
    """

rule annotate_variants:
  input:
    calls="results/{sample}/variants/{vartype}.vcf.gz",
    idx="results/{sample}/variants/{vartype}.vcf.gz.tbi",
    fasta="resources/refs/genome.fasta",
    cache="resources/vep/cache",
    plugins="resources/vep/plugins",
    wt_plugin="resources/vep/plugins/Wildtype.pm"
  output:
    calls="results/{sample}/annotation/{vartype}.vcf",
    stats="results/{sample}/annotation/{vartype}.html"
  params:
    plugins=["NMD","Wildtype","Downstream"],
    extra="--everything",  # optional: extra arguments
  log:
    "logs/vep/{sample}_{vartype}_annotate.log",
  threads: 4
  wrapper:
    "v1.31.1/bio/vep/annotate"


