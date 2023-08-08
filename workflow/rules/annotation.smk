rule download_vep_plugins:
    output:
        directory("resources/vep/plugins")
    params:
        release=110
    wrapper:
        "v1.31.1/bio/vep/plugins"

rule download_wildtype_plugin:
  output:
    "resources/vep/plugins/Wildtype.pm"
  conda:
    "../envs/basic.yml"
  log:
    "logs/vep/plugins.log"
  shell:
    "curl -L https://raw.githubusercontent.com/griffithlab/pVACtools/v4.0.1/pvactools/tools/pvacseq/VEP_plugins/Wildtype.pm > resources/vep/plugins/Wildtype.pm"



#rule download_vep_cache:
  #output:
    #directory("resources/vep/cache")
  #conda:
    #"../envs/basic.yaml"
  #log:
    #"logs/vep/cache.log"
  #shell:
    #"""
      #wget https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
      #tar xzf homo_sapiens_vep_110_GRCh38.tar.gz
      #mv homo_sapiens_vep_110_GRCh38 resources/vep/cache
    #"""

#rule get_vep_cache:
  #output:
      #directory("resources/vep/cache"),
  #params:
      #species="homo_sapiens",
      #build="GRCh38",
      #release="110",
  #log:
      #"logs/vep/cache.log",
  #cache: "omit-software"  # save space and time with between workflow caching (see docs)
  #wrapper:
      #"v2.2.1/bio/vep/cache"


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


rule variants_to_peptides:
  input:
    get_variants,
  output:
    "results/{sample}/neoantigens/result.tsv"
  log:
    "logs/vep/{sample}_variants_to_peptides.log"
  conda:
    "../envs/variants_to_peptides.yml"
  shell:
    """
      python3 workflow/scripts/variants_to_peptide.py \
          resources/refs/Homo_sapiens/GRCh38.pep.all.fa \
          {input} {output}
    """






















# annotation of long indels (as determined by transindel)
#rule vep_longindels:
  #input:
    #calls="results/{sample}/variants/long.indels.vcf.gz",
    #idx="results/{sample}/variants/long.indels.vcf.gz.tbi",
    #cache="resources/vep/cache",
    #plugins="resources/vep/plugins"
  #output:
    #calls="results/{sample}/annotation/long.indels.vcf",
    #stats="results/{sample}/annotation/long.indels.html"
  #params:
    #plugins=["NMD","Wildtype","Downstream"],
    #extra="--everything",  # optional: extra arguments
  #log:
    #"logs/vep/{sample}_annotate_long_indel.log",
  #threads: 4
  #wrapper:
    #"v2.2.1/bio/vep/annotate"

# annotation of short indels (as determined by GATK)
#rule vep_somatic_short_indel:
  #input:
    #calls="results/{sample}/variants/somatic.short.indels.vcf",
    #plugins="resources/vep/plugins",
    #fasta="resources/refs/genome.fasta",
    #fai="resources/refs/genome.fasta.fai", # fasta index
    #gff="resources/refs/genome.gtf.gz",
    #csi="resources/refs/genome.gtf.gz.csi",
  #output:
    #calls="results/{sample}/annotation/somatic.short.indels.vcf",
    #stats="results/{sample}/annotation/somatic.short.indels.html"
  #params:
    #plugins=["NMD"],
    #extra="--everything",  # optional: extra arguments
  #log:
    #"logs/vep/{sample}_annotate.log",
  #threads: 4
  #wrapper:
    #"v1.31.1/bio/vep/annotate"

## annotation of snvs (as determined by GATK)
#rule vep_somatic_snvs:
  #input:
    #calls="results/{sample}/variants/somatic.snvs.vcf",
    #plugins="resources/vep/plugins",
    #fasta="resources/refs/genome.fasta",
    #fai="resources/refs/genome.fasta.fai", # fasta index
    #gff="resources/refs/genome.gtf.gz",
    #csi="resources/refs/genome.gtf.gz.csi",
  #output:
    #calls="results/{sample}/annotation/somatic.snvs.vcf",
    #stats="results/{sample}/annotation/somatic.snvs.html"
  #params:
    #plugins=["NMD"],
    #extra="--everything",  # optional: extra arguments
  #log:
    #"logs/vep/{sample}_annotate.log",
  #threads: 4
  #wrapper:
    #"v1.31.1/bio/vep/annotate"










