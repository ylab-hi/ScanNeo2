rule spladder:
    input:
        bam = "results/{sample}/rnaseq/align/ready.bam",
        bamidx = "results/{sample}/rnaseq/align/ready.bam.bai"
    output:
        directory("results/{sample}/rnaseq/altsplicing/spladder/")
    conda: 
        "../envs/spladder.yml"
    log:
        "logs/{sample}/spladder/build.log"
    params:
    shell:
        """
        spladder build -b {input.bam} \
        -a {config[annotation]} -o {output} \
        --filter-overlap-exons --no-primary-only \
        --confidence 2 --quantify-graph \
        --qmode all > {log} 2>&1
        """

rule splicing_to_vcf:
    input:
        "results/as/{sample}/merge_graphs.txt"
    output:
        "results/as/{sample}.vcf"
    conda:
        "../envs/spladder.yml"
    shell:
        """
        """

rule combine_splicing:
    input:
        expand("results/as/{sample}.vcf", sample=config["rnaseq"])
    output:
        "results/as/all.vcf"
    shell:
        """
        """
