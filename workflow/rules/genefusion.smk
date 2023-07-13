rule arriba:
    input:
        bam="results/{sample}/rnaseq/align/ready.bam",
        genome="resources/refs/genome.fasta",
        annotation="resources/refs/genome.gtf",
        custom_blacklist=[],
    output:
        fusions="results/{sample}/genefusion/fusions.tsv",
        discarded="results/{sample}/genefusion/fusions.discarded.tsv",  # optional
    log:
        "logs/{sample}/arriba/run.log",
    params:
        genome_build="GRCh38",
        default_blacklist=False,  # optional
        default_known_fusions=True,  # optional
        sv_file="",  # optional
        extra="-i 1,2 -G 'gene_name=gene gene_id=gene_id transcript_id=transcript_id feature_exon=exon feature_CDS=CDS'"
    threads: 8
    wrapper:
        "v1.29.0/bio/arriba"

rule fusions_to_vcf:
    input:
        "results/{sample}/genefusion/fusions.tsv",
    output:
        "results/{sample}/variants/fusions.vcf"
    log:
        "logs/fusion_to_vcf_{sample}.log"
    shell:
        "workflow/scripts/convert_fusions_to_vcf.sh resources/refs/genome.fasta {input} {output} > {log} 2>&1"

