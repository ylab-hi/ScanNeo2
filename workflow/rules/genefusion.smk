rule arriba:
    input:
        bam = "results/{sample}/rnaseq/align/{group}_final_STAR.bam",
        genome="resources/refs/genome.fasta",
        annotation="resources/refs/genome.gtf",
        custom_blacklist=[],
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
        extra=f"""-G 'gene_name=gene gene_id=gene_id transcript_id=transcript_id feature_exon=exon feature_CDS=CDS' \
          -E {config['genefusion']['maxevalue']} \
          -S {config['genefusion']['suppreads']} \
          -L {config['genefusion']['maxidentity']} \
          -H {config['genefusion']['hpolymerlen']} \
          -R {config['genefusion']['readthroughdist']} \
          -A {config['genefusion']['minanchorlen']} \
          -M {config['genefusion']['splicedevents']} \
          -K {config['genefusion']['maxkmer']} \
          -F {config['genefusion']['fraglen']} \
          -V {config['genefusion']['maxmismatch']} \
          -U {config['genefusion']['maxsuppreads']}"""
    threads: config['threads']
    wrapper:
        "v1.29.0/bio/arriba"

rule fusions_to_vcf:
    input:
        "results/{sample}/rnaseq/genefusion/{group}_fusions.tsv",
    output:
        "results/{sample}/rnaseq/genefusion/{group}_fusions.vcf"
    log:
        "logs/{sample}/genefusion/fusion_to_vcf_{group}.log"
    shell:
        "workflow/scripts/convert_fusions_to_vcf.sh resources/refs/genome.fasta {input} {output} > {log} 2>&1"

rule combine_fusions:
  input:
    expand("results/{sample}/rnaseq/genefusion/{group}_fusions.vcf", 
           sample=config['data']['name'],
           group=list(config['data']['rnaseq'].keys()))
  output:
    "results/{sample}/variants/fusions.vcf"
  message:
    "Combining somatic SNVs/Indels with Mutect2 on sample:{wildcards.sample}"
  log: 
    "logs/{sample}/combine/fusions.vcf"
  conda:
    "../envs/manipulate_vcf.yml"
  shell:
    """
      python workflow/scripts/combine_vcf.py '{input}' {output} > {log} 2>&1
    """