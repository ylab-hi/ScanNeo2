rule arriba:
    input:
        bam="results/{sample}/rnaseq/align/{group}_final_STAR.bam",
        genome="resources/refs/genome.fasta",
        annotation="resources/refs/genome.gtf",
    output:
        fusions="results/{sample}/rnaseq/genefusion/{group}_fusions.tsv",
        discarded="results/{sample}/rnaseq/genefusion/{group}_fusions.discarded.tsv",
    log:
        "logs/{sample}/genefusion/arriba_{group}.log",
    threads: config["threads"]
    params:
        genome_build="GRCh38",
        # arriba's blacklist is its main artifact filter; without it arriba
        # retains a huge candidate set (tens of millions of alignments) and
        # OOMs. The GRCh38 blacklist ships with the arriba conda package.
        default_blacklist=True,
        default_known_fusions=True,
        sv_file="",
        extra=_arriba_extra(),
    wrapper:
        "v2.6.0/bio/arriba"
