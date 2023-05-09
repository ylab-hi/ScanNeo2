rule scanexitron:
    input: 
        "results/preproc/post/{sample}_dedup.bam",
    output:
        "results/exitron/{sample}/res.exitron",
    log:
        "logs/exitron_{sample}.log"
    conda:
        "../envs/scanexitron.yml",
    shell:
        "python3 workflow/scripts/ScanExitron/ScanExitron.py -i {input} -r hg38 > {output}"
        
rule exitron_to_vcf:
    input:
        "results/exitron/{sample}/res.exitron"
    output:
        "results/exitron/{sample}/res.vcf"
    log:
        "logs/exitron2vcf_{sample}.log"
    shell:
        "workflow/scripts/exitron2vcf.py {config[refgen]} {input} {output} 2> logs/error.err"

rule combine_exitrons:
    input:
        expand("results/exitron/{sample}/res.vcf", sample=config["rnaseq"])
    output:
        "results/exitron/all.vcf"
    shell:
        """
        """
