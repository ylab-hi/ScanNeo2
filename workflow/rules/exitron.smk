rule prepare_scanexitron_cds:
    input:
        gtf="resources/refs/genome.gtf",
    output:
        # ScanExitron derives this path from the GTF (<gtf_stem>.CDS.bed) and
        # caches it there. Its cache writer is not concurrency-safe (bare
        # exists-check + in-place write), so pre-generate it once here and have
        # every scanexitron job depend on it -- concurrent samples then reuse
        # the cache instead of racing to regenerate it.
        "resources/refs/genome.CDS.bed",
    log:
        "logs/ref/prepare_scanexitron_cds.log",
    conda:
        "../envs/scanexitron.yml"
    shell:
        """
        PYTHONPATH=workflow/scripts/scanexitron/src python3 -c \
            'import sys; from scanexitron.gtf import extract_cds_bed; extract_cds_bed(sys.argv[1])' \
            {input.gtf} >{log} 2>&1
        """


rule scanexitron:
    input:
        bam="results/{sample}/rnaseq/align/{group}_final_STAR.bam",
        idx="results/{sample}/rnaseq/align/{group}_final_STAR.bam.bai",
        fasta="resources/refs/genome.fasta",
        gtf="resources/refs/genome.gtf",
        cds="resources/refs/genome.CDS.bed",
    output:
        "results/{sample}/rnaseq/exitron/{group}.exitron",
    log:
        "logs/{sample}/exitron/scanexitron_{group}.log",
    conda:
        "../envs/scanexitron.yml"
    threads: config["threads"]
    params:
        mapq=config["mapq"],
        ao=config["exitronsplicing"]["ao"],
        pso=config["exitronsplicing"]["pso"],
        strand=config["exitronsplicing"]["strand"],
    message:
        "Detect exitrons on sample:{wildcards.sample} of group:{wildcards.group}"
    shell:
        # ScanExitron 1.4.0 typer CLI (vendored submodule); it writes
        # {output-prefix}.exitron and keeps intermediates in a TemporaryDirectory
        # (the fix branch), so no CWD cleanup is needed.
        """
        PYTHONPATH=workflow/scripts/scanexitron/src python3 -m scanexitron run \
            -i {input.bam} \
            -r {input.fasta} \
            -g {input.gtf} \
            -a {params.ao} \
            -p {params.pso} \
            -m {params.mapq} \
            -s {params.strand} \
            -t {threads} \
            -o results/{wildcards.sample}/rnaseq/exitron/{wildcards.group} >{log} 2>&1
        """


rule exitron_to_vcf:
    input:
        "results/{sample}/rnaseq/exitron/{group}.exitron",
    output:
        "results/{sample}/rnaseq/exitron/{group}_exitrons.vcf",
    log:
        "logs/{sample}/exitron/exitron_to_vcf_{group}.log",
    conda:
        "../envs/manipulate_vcf.yml"
    shell:
        """
        python workflow/scripts/exitron2vcf.py \
            {input} {output} \
            resources/refs/genome.fasta >{log} 2>&1
        """


rule exitron_augment:
    input:
        "results/{sample}/rnaseq/exitron/{group}_exitrons.vcf",
    output:
        "results/{sample}/rnaseq/exitron/{group}_exitrons_augmented.vcf",
    log:
        "logs/{sample}/exitron/exitron_augment_{group}.log",
    conda:
        "../envs/manipulate_vcf.yml"
    message:
        "Augmenting exitrons on sample:{wildcards.sample} of group:{wildcards.group}"
    shell:
        """
        python workflow/scripts/add_infos_to_vcf.py \
            {input} \
            exitron \
            {wildcards.group} \
            {output} >{log} 2>&1
        """


rule sort_exitron:
    input:
        "results/{sample}/rnaseq/exitron/{group}_exitrons_augmented.vcf",
    output:
        "results/{sample}/rnaseq/exitron/{group}_exitrons.vcf.gz",
    log:
        "logs/{sample}/exitron/sort_exitron_{group}.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Sorting and compressing exitrons on sample:{wildcards.sample} of group:{wildcards.group}"
    shell:
        """
        (bcftools sort {input} -o - | bcftools view -O z -o {output}) >{log} 2>&1
        """


rule combine_exitrons:
    input:
        get_exitrons,
    output:
        "results/{sample}/variants/exitrons.vcf.gz",
    log:
        "logs/{sample}/exitron/combine_exitrons.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Combining exitrons on sample:{wildcards.sample}"
    shell:
        """
        (bcftools concat --naive-force -O z {input} -o - | bcftools sort -O z -o {output}) >{log} 2>&1
        """
