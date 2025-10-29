rule download_vep_plugins:
  output:
    directory("resources/vep/plugins")
  message:
    "Downloading VEP plugins"
  log: 
    "logs/vep/download_plugins.log"
  conda:
    "../envs/basic.yml"
  params:
    release=f"""{config['reference']['release']}"""
  shell:
    """
      mkdir -p resources/vep/plugins/
      curl -L -o {output}/NMD.pm https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/{params.release}/NMD.pm
      curl -L -o {output}/Downstream.pm https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/{params.release}/Downstream.pm
      curl -L -o {output}/Wildtype.pm https://raw.githubusercontent.com/griffithlab/pVAC-Seq/master/pvacseq/VEP_plugins/Wildtype.pm
    """

rule download_vep_cache:
  output:
    directory("resources/vep/cache")
  message:
    "Downloading VEP cache"
  log:
    "logs/vep/cache.log"
  conda:
    "../envs/basic.yml"
  params:
    release=f"""{config['reference']['release']}"""
  shell:
    """
      mkdir -p {output}
      curl -L -C - \
          --retry 10 \
          --retry-delay 10 \
          --retry-all-errors \
          -o resources/vep/cache/homo_sapiens_vep_{params.release}_GRCh38.tar.gz \
          https://ftp.ensembl.org/pub/release-{params.release}/variation/indexed_vep_cache/homo_sapiens_vep_{params.release}_GRCh38.tar.gz

      tar -xzf resources/vep/cache/homo_sapiens_vep_{params.release}_GRCh38.tar.gz -C resources/vep/cache/
      rm resources/vep/cache/homo_sapiens_vep_{params.release}_GRCh38.tar.gz

    """
#      curl -L https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz \
          # https://g-a8b222.dd271.03c0.data.globus.org/ensemblorg/pub/release-{params.release}/variation/indexed_vep_cache/homo_sapiens_vep_{params.release}_GRCh38.tar.gz

    

rule index_variants:
  input:
    "results/{sample}/variants/{vartype}.vcf.gz"
  output:
    "results/{sample}/variants/{vartype}.vcf.gz.tbi"
  log:
    "logs/indexvcf/{sample}_{vartype}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      tabix {input}
    """

rule annotate_variants:
  input:
    calls="results/{sample}/variants/{vartype}.vcf.gz",
    idx="results/{sample}/variants/{vartype}.vcf.gz.tbi",
    fasta="resources/refs/genome.fasta",
    cache="resources/vep/cache",
    plugins="resources/vep/plugins"
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
    "v5.9.0/bio/vep/annotate"


