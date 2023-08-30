rule get_genome:
  output:
    genome="resources/refs/genome.fasta",
    annotation="resources/refs/genome.gtf"
  message:
    "Download reference genome and annotation"
  conda:
    "../envs/basic.yml"
  log:
    "logs/get-genome.log",
  params:
  shell:
    """
      mkdir -p resources/refs
      curl -L -o {output.genome}.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz 
      gzip -d {output.genome}.gz
      curl -L -o {output.annotation}.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz 
      gzip -d {output.annotation}
    """

      #curl -L -o - https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
        #| gzip -d - > {output.genome}
      #curl -L -o - https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz \
          #| gzip -d - > {output.annotation}

rule genome_index: 
  input:
    "resources/refs/genome.fasta"
  output:
    "resources/refs/genome.fasta.fai",
  message:
    "Create index of reference genome"
  log:
    "logs/samtools/index_genome.log",
  params:
    extra="",  # optional params string
  wrapper:
    "v1.28.0/bio/samtools/faidx"

rule annotation_sort_bgzip:
  input:
    "resources/refs/genome.gtf"
  output:
    "resources/refs/genome.gtf.gz"
  message:
    "Sort of bgzip the annotations file"
  log:
    "logs/bgzip_annotation"
  conda:
    "../envs/basic.yml"
  shell:
    """
      (grep "^#" {input}; grep -v "^#" {input} | sort -t"`printf '\t'`" -k1,1 -k4,4n) | bgzip > {output}
    """

rule tabix:
  input:
      "resources/refs/genome.gtf.gz",
  output:
      "resources/refs/genome.gtf.gz.csi",
  log:
      "logs/tabix/annotation.log",
  params:
      # pass arguments to tabix (e.g. index a vcf)
      "-C -p gff",
  wrapper:
      "v1.29.0/bio/tabix/index"

rule star_index:
  input:
      fasta="resources/refs/genome.fasta"
  output:
      directory("resources/refs/star/"),
  message:
      "Creating STAR index"
  threads: config['threads']
  params:
      extra="--genomeSAindexNbases 5 ",
  log:
      "logs/star/create_star_index.log",
  wrapper:
      "v1.26.0/bio/star/index"

rule bwa_index:
  input:
    fasta="resources/refs/genome.fasta"
  output:
    idx = multiext("resources/refs/bwa/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
  log:
    "logs/bwa_index.log",
  params:
    algorithm="bwtsw",
  wrapper:
    "v1.26.0/bio/bwa/index"


rule create_sequence_dictionary:
  input:
    "resources/refs/genome.fasta"
  output:
    "resources/refs/genome.dict"
  message:
    "Create sequence dictionary of reference genome"
  log:
    "logs/picard/create_dict.log"
  params:
    extra="",  # optional: extra arguments for picard.
  resources:
    mem_mb=10024,
  wrapper:
    "v1.31.1/bio/picard/createsequencedictionary"
