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
    shell:
        """
          spladder build -b {input.bam} \
              -a resources/refs/genome.gtf \
              -o {output} --filter-overlap-exons \
              --no-primary-only --confidence {config[altsplicing][confidence]} \
              --iterations {config[altsplicing][iterations]} \
              --ase-edge-limit {config[altsplicing][edgelimit]} \
              --quantify-graph --qmode all > {log} 2>&1
        """

rule splicing_to_vcf:
    input:
      "results/{sample}/rnaseq/altsplicing/spladder/{group}/merge_graphs_{type}_{confidence}.confirmed.txt.gz"
    output:
        "results/as/{sample}.vcf"
    conda:
        "../envs/spladder.yml"
    shell:
        """
        """

#rule combine_splicing:
    #input:
        #expand("results/as/{sample}.vcf", sample=config["rnaseq"])
    #output:
        #"results/as/all.vcf"
    #shell:
        #"""
        #"""
