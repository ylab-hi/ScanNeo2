rule download_vep_plugins:
    output:
        directory("resources/vep/plugins")
    params:
        release=110
    wrapper:
        "v1.31.1/bio/vep/plugins"

rule download_vep_cache:
  output:
    directory("resources/vep/cache")
  conda:
    "../envs/basic.yml"
  log:
    "logs/vep/cache.log"
  shell:
    """
      curl -L -o https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz 
      | tar xvf - -C resources/vep/cache
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
    cache="resources/vep/cache",
    plugins="resources/vep/plugins",
  output:
    calls="results/{sample}/annotation/{vartype}.vcf",
    stats="results/{sample}/annotation/{vartype}.html"
  params:
    plugins=["NMD","Wildtype","Downstream"],
    extra="--everything --plugin Downstream --plugin Wildtype --plugin NMD",  # optional: extra arguments
  log:
    "logs/vep/{sample}_{vartype}_annotate.log",
  threads: 4
  wrapper:
    "v1.31.1/bio/vep/annotate"


