rule spladder:
    input:
        bam="results/{sample}/rnaseq/align/{group}_final_STAR.bam",
        bamidx="results/{sample}/rnaseq/align/{group}_final_STAR.bam.bai",
        gtf="resources/refs/genome.gtf",
    output:
        directory("results/{sample}/rnaseq/altsplicing/spladder_run/{group}"),
    log:
        "logs/{sample}/altsplicing/spladder_{group}.log",
    conda:
        "../envs/spladder.yml"
    threads: config["threads"]
    params:
        confidence=f"""{config["altsplicing"]["confidence"]}""",
        iteration=f"""{config["altsplicing"]["iterations"]}""",
        edgelimit=f"""{config["altsplicing"]["edgelimit"]}""",
    shell:
        """
        bash workflow/scripts/run_spladder.sh \
            {input.bam} \
            {threads} \
            {input.gtf} \
            {output} \
            {params.confidence} \
            {params.iteration} \
            {params.edgelimit} \
            {log}
        """
