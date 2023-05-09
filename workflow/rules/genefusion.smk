rule arriba:
    input:
        bam="results/preproc/post/{sample}_dedup.bam",
        genome="refs/genome.fasta",
        annotation="refs/genome.gtf",
        custom_blacklist=[],
    output:
        fusions="results/genefusion/{sample}/fusions.tsv",
        discarded="results/genefusion/{sample}/fusions.discarded.tsv",  # optional
    log:
        "logs/arriba_{sample}.log",
    params:
        genome_build="GRCh38",
        default_blacklist=False,  # optional
        default_known_fusions=True,  # optional
        sv_file="",  # optional
        extra="-i 1,2 -G 'gene_name=gene_name gene_id=gene_id transcript_id=transcript_id feature_exon=exon feature_CDS=CDS'"
    threads: 8
    wrapper:
        "v1.29.0/bio/arriba"

rule fusions_to_vcf:
    input:
        "results/genefusion/{sample}/fusions.tsv"
    output:
        "results/genefusion/{sample}/fusions.vcf"
    log:
        "logs/fusion_to_vcf_{sample}.log"
    shell:
        "workflow/scripts/convert_fusions_to_vcf.sh {config[refgen]} {input} {output} 2> {log}"


rule combine_fusions:
    input:
        expand( "results/genefusion/{sample}/fusions.vcf", sample=config["rnaseq"])
    output:
        "results/genefusion/all.vcf"
    shell:
        """
        """
