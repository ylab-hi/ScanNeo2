rule get_genome:
    output:
        genome="refs/genome.fasta",
        annotation="refs/genome.gtf"
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    shell:
        """
        curl -L https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz | gzip -d > {output.genome}
        curl -L https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz | gzip -d > {output.annotation}
        """

rule genome_index: 
    input:
        "refs/genome.fasta"
    output:
        "refs/genome.fasta.fai",
    log:
        "logs/samtools/index_genome.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.28.0/bio/samtools/faidx"


rule annotation_sort_bgzip:
    input:
        "refs/genome.gtf"
    output:
        "refs/genome.gtf.gz"
    shell:
        """
            (grep "^#" {input}; grep -v "^#" {input} | sort -t"`printf '\t'`" -k1,1 -k4,4n) | bgzip > {output}
        """

rule tabix:
    input:
        "refs/genome.gtf.gz",
    output:
        "refs/genome.gtf.gz.csi",
    log:
        "logs/tabix/annotation.log",
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-C -p gff",
    wrapper:
        "v1.29.0/bio/tabix/index"

rule star_index:
    input:
        "refs/genome.fasta"
    output:
        directory("refs/star/"),
    message:
        "Creating STAR index",
    threads: 60
    params:
        extra="--genomeSAindexNbases 5 ",
    log:
        "logs/star/create_star_index.log",
    wrapper:
        "v1.26.0/bio/star/index"
