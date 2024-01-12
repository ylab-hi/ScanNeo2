rule arriba:
    input:
        bam = "results/{sample}/rnaseq/align/{group}_final_STAR.bam",
        genome="resources/refs/genome.fasta",
        annotation="resources/refs/genome.gtf",
#        custom_blacklist=[],
    output:
        fusions="results/{sample}/rnaseq/genefusion/{group}_fusions.tsv",
        discarded="results/{sample}/rnaseq/genefusion/{group}_fusions.discarded.tsv",  # optional
    log:
        "logs/{sample}/arriba/{group}.log",
    params:
        genome_build="GRCh38",
        default_blacklist=False,  # optional
        default_known_fusions=True,  # optional
        sv_file="",  # optional
        extra=f"""-G 'gene_name=gene_name|gene_id gene_id=gene_id transcript_id=transcript_id feature_exon=exon feature_CDS=CDS' \
          -E {config['genefusion']['maxevalue']} \
          -S {config['genefusion']['suppreads']} \
          -L {config['genefusion']['maxidentity']} \
          -H {config['genefusion']['hpolymerlen']} \
          -R {config['genefusion']['readthroughdist']} \
          -A {config['genefusion']['minanchorlen']} \
          -M {config['genefusion']['splicedevents']} \
          -K {config['genefusion']['maxkmer']} \
          -F {config['genefusion']['fraglen']} \
          -V {config['genefusion']['maxmismatch']}
        """
    threads: config['threads']
    wrapper:
        "v2.6.0/bio/arriba"
          
