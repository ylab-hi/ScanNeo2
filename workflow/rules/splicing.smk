def comma_join(file_list):
    return ",".join(file_list)

rule spladder:
    input:
        bam = expand("results/preproc/post/{sample}_dedup.bam", sample=config["samples"]),
        bamidx = expand("results/preproc/post/{sample}_dedup.bam.bai", sample=config["samples"]),
    output:
       directory("results/splicing/spladder") 
    conda: 
        "../envs/spladder.yml"
    log:
        "logs/spladder.log"
    params:
        bams = comma_join(expand("results/preproc/post/{sample}_dedup.bam", sample=config["samples"])),
    shell:
        """
        spladder build -b {params.bams} \
        -a {config[annotation]} -o {output} \
        --filter-overlap-exons --no-primary-only \
        --confidence 2 --quantify-graph \
        --qmode all 2> {log}
        """

