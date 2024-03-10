rule get_genome:
  output:
    genome="resources/refs/genome_tmp.fasta",
    annotation="resources/refs/genome_tmp.gtf",
    peptide="resources/refs/peptide.fasta"
  message:
    "Download reference genome and annotation"
  conda:
    "../envs/basic.yml"
  log:
    "logs/get-genome.log",
  params:
    release=f"""{config['reference']['release']}"""
  shell:
    """
      curl -L -o {output.genome}.gz https://ftp.ensembl.org/pub/release-{params.release}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 
      gzip -d {output.genome}.gz
      curl -L -o {output.annotation}.gz https://ftp.ensembl.org/pub/release-{params.release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz 
      gzip -d {output.annotation}.gz
      curl -L -o {output.peptide}.gz https://ftp.ensembl.org/pub/release-{params.release}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
      gzip -d {output.peptide}.gz
    """

rule mod_genome:
  input:
    genome="resources/refs/genome_tmp.fasta",
    annotation="resources/refs/genome_tmp.gtf"
  output:
    genome="resources/refs/genome.fasta",
    annotation="resources/refs/genome.gtf"
  message:
    "Modify header in reference genome/annotation"
  conda:
    "../envs/basic.yml"
  log:
    "logs/mod-genome.log"
  shell:
    """
      python3 workflow/scripts/reference/modify_ensembl_header.py {input.genome} {output.genome} {input.annotation} {output.annotation}
    """

# this forces to redownload the reference on each execution
##      rm resources/refs/genome_tmp.fasta
#      rm resources/refs/genome_tmp.gtf

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


rule get_hla_info:
  output:
    "resources/hla/hla_gen.fasta"
  message:
    "Download full resolution HLA information"
  log:
    "logs/get_hla_info.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      curl -L -o {output} ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta
    """

rule create_hla_idx_bowtie:
  input:
    ref="resources/hla/hla_gen.fasta"
  output:
    multiext(
            "resources/hla/bowtie2_index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
  message:
    "Create bowtie index of HLA reference"
  log:
    "logs/bowtie2/create_hla_idx.log"
  params:
      extra="",  # optional parameters
  threads: 8
  wrapper:
      "v2.11.1/bio/bowtie2/build"
