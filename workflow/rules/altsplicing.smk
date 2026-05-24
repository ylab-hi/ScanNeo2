rule spladder:
    input:
        bam="results/{sample}/rnaseq/align/{group}_final_STAR.bam",
        bamidx="results/{sample}/rnaseq/align/{group}_final_STAR.bam.bai",
    output:
        directory("results/{sample}/rnaseq/altsplicing/spladder/{group}"),
    log:
        "logs/{sample}/altsplicing/spladder_{group}.log",
    conda:
        "../envs/spladder.yml"
    threads: 20
    params:
        confidence=f"""{config["altsplicing"]["confidence"]}""",
        iteration=f"""{config["altsplicing"]["iterations"]}""",
        edgelimit=f"""{config["altsplicing"]["edgelimit"]}""",
    shell:
        """
        bash workflow/scripts/run_spladder.sh \
            {input.bam} \
            {threads} \
            resources/refs/genome.gtf \
            {output} \
            {params.confidence} \
            {params.iteration} \
            {params.edgelimit} \
            {log}
        """


rule splicing_to_vcf:
    input:
        "results/{sample}/rnaseq/altsplicing/spladder/{group}",
    output:
        "results/{sample}/rnaseq/altsplicing/spladder/{group}_altsplicing.vcf",
    log:
        "logs/{sample}/altsplicing/splicing_to_vcf_{group}.log",
    conda:
        "../envs/manipulate_vcf.yml"
    message:
        "Converting splicing events to VCF format"
    shell:
        """
        python workflow/scripts/altsplc2vcf.py \
            -i {input} -r resources/refs/genome.fasta \
            -g {wildcards.group} -o {output} >{log} 2>&1
        """


rule sort_altsplicing:
    input:
        "results/{sample}/rnaseq/altsplicing/spladder/{group}_altsplicing.vcf",
    output:
        "results/{sample}/rnaseq/altsplicing/spladder/{group}_altsplicing.vcf.gz",
    log:
        "logs/{sample}/altsplicing/sort_altsplicing_{group}.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Sorting and compressing splicing events on sample:{wildcards.sample} of group:{wildcards.group}"
    shell:
        """
        (bcftools sort {input} -o - | bcftools view -O z -o {output}) >{log} 2>&1
        """


rule combine_altsplicing:
    input:
        get_altsplicing,
    output:
        "results/{sample}/variants/altsplicing.vcf.gz",
    log:
        "logs/{sample}/altsplicing/combine_altsplicing.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Combining alternative splicing events on sample:{wildcards.sample}"
    shell:
        """
        (bcftools concat --naive-force -O z {input} -o - | bcftools sort -O z -o {output}) >{log} 2>&1
        """
