rule spladder:
    input:
        bam = "results/{sample}/rnaseq/align/{group}_final_STAR.bam",
        bamidx = "results/{sample}/rnaseq/align/{group}_final_STAR.bam.bai"
    output:
        directory("results/{sample}/rnaseq/altsplicing/spladder/{group}")
    conda: 
        "../envs/spladder.yml"
    log:
        "logs/{sample}/spladder/{group}_build.log"
    params:
      confidence = "--confidence {config[altsplicing][confidence]}",
      iteration = "--iterations {config[altsplicing][iterations]}",
      edgelimit = "--ase-edge-limit {config[altsplicing][edgelimit]}",
    shell:
        """
          spladder build -b {input.bam} \
              -a resources/refs/genome.gtf \
              -o {output} --filter-overlap-exons \
              --no-primary-only --quantify-graph \
              --qmode all > {log} 2>&1
        """

#rule splicing_to_vcf:
    #input:
      #"results/{sample}/rnaseq/altsplicing/spladder/{group}/merge_graphs_{type}_{confidence}.confirmed.txt.gz"
    #output:
        #"results/as/{sample}.vcf"
    #message:
      #"Converting splicing events to VCF format"
    #log:
      #"logs/{sample}/spladder/{group}_vcf.log"
    #conda:
        #"../envs/spladder.yml"
    #shell:
        #"""
        #"""

#rule combine_splicing:
    #input:
        #expand("results/as/{sample}.vcf", sample=config["rnaseq"])
    #output:
        #"results/as/all.vcf"
    #shell:
        #"""
        #"""
