rule get_hla_panel:
    output:
        dna="misc/hla/hla_reference_dna.fasta",
        rna="misc/hla/hla_reference_rna.fasta"
    shell:
        """
        curl -o {output.dna} https://raw.githubusercontent.com/FRED-2/OptiType/v1.3.5/data/hla_reference_dna.fasta
        curl -o {output.rna} https://raw.githubusercontent.com/FRED-2/OptiType/v1.3.5/data/hla_reference_rna.fasta
        """

rule index_hla_panel:
    input:
        dna="misc/hla/hla_reference_dna.fasta",
        rna="misc/hla/hla_reference_rna.fasta",
    output:
        "misc/hla/hla_dna.index.lf.drp", 
        "misc/hla/hla_rna.index.lf.drp", 
    log:
        "logs/yara_indexer.log"
    conda:
        "../envs/yara.yml"
    shell:
        """
        yara_indexer -o misc/hla/hla.index {input.dna} > {log}
        yara_indexer -o misc/hla/hla.index {input.rna} >> {log}
        """

rule prepare_bams:
    input:
        expand("results/preproc/align/{sample}.bam", sample=config["rnaseq"])
    output:
        "results/hla/all.fastq"
    log:
        "logs/bam2fastq.log"
    conda:
        "../envs/yara.yml"
    shell:
        "samtools merge -f - {input} | samtools fastq - > {output}"

rule filter_hla:
    input:
        reads="results/hla/all.fastq",
        yidx="misc/hla/hla.index.lf.drp"
    output:
        "results/hla/all_hla.fastq"
    conda:
        "../envs/yara.yml"
    shell:
        """
        yara_mapper -e 3 -f bam -u misc/hla/hla.index {input.reads} \
        | samtools view -F 4 -h -b1 - | samtools bam2fq - > {output}
        """

rule hla_genotyping:
    input:
        reads = "results/hla/all_hla.fastq"
    output:
        pdf="results/hla/all_coverage_plot.pdf",
        tsv="results/hla/all_result.tsv",
    conda: 
        "../envs/optitype.yml"
    log:
        "logs/optitype/call.log"
    params:
        # Type of sequencing data. Can be 'dna' or 'rna'. Default is 'dna'.
        sequencing_type="rna",
        # optiype config file, optional
        config="",
        # additional parameters
        extra=""
    wrapper:
        "v1.26.0/bio/optitype"
        
