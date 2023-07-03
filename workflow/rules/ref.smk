rule get_genome:
    output:
        genome="resources/refs/genome.fasta",
        annotation="resources/refs/genome.gtf"
    conda:
        "../envs/basic.yml"
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    shell:
        """
        curl -L -o - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz | gzip -d - >{output.genome}
        curl -L -o - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz | gzip -d - > {output.annotation}       
        """

rule genome_index: 
    input:
        "resources/refs/genome.fasta"
    output:
        "resources/refs/genome.fasta.fai",
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
        "Creating STAR index",
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
