rule download_vep_plugins:
  output:
    directory("resources/vep/plugins")
  message:
    "Downloading VEP plugins"
  log: 
    "logs/vep/download_plugins.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      mkdir -p resources/vep/plugins/
      curl -L -o {output}/NMD.pm https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/NMD.pm
      curl -L -o {output}/Downstream.pm https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/110/Downstream.pm
      curl -L -o {output}/Wildtype.pm https://raw.githubusercontent.com/griffithlab/pVAC-Seq/master/pvacseq/VEP_plugins/Wildtype.pm
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
      curl -L https://g-a8b222.dd271.03c0.data.globus.org/ensemblorg/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz \
      | tar -xz -C resources/vep/cache
    """
#      curl -L https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz \

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
    "v1.31.1/bio/vep/annotate"


