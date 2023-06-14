# if the case of specified bam files, convert to fastq
if rnaseq_filetype == ".bam":
    rule bam_to_fastq:
        input:
            get_bams,
        output:
            "results/{sample}/rnaseq/preproc/pre/inputreads.fq.gz"
        conda:
            "../envs/samtools.yml"
        log:
            "logs/samtools/bam2fastq/{sample}.log"
        threads: config['threads']
        shell:
            """
                samtools collate -Oun128 -@ {threads} {input} \
                | samtools fastq -OT RG,BC -@ {threads} - | gzip -c - > {output}
            """


checkpoint split_fastq:
    input:
        "results/{sample}/rnaseq/preproc/pre/inputreads.fq.gz"
    output:
        directory("results/{sample}/rnaseq/preproc/pre/fqs/")
    log:
        "logs/split_fastq/{sample}.log"
    conda:
        "../envs/splitfastq.yml"
    shell:
        """
            python workflow/scripts/splitfastq.py {input} \
            results/{wildcards.sample}/rnaseq/preproc/pre/ \
            10000000 > {log} 2>&1
        """

#extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --limitBAMsortRAM 11209980843"
rule star_align:
    input:
        fq1 = "results/{sample}/rnaseq/preproc/pre/fqs/inputreads_{i}.fq.gz",
        idx = "resources/refs/star/",
    output:
        aln = "results/{sample}/rnaseq/preproc/align/bams/aln_{i}.bam",
        log = "results/{sample}/rnaseq/preproc/align/bams/star.align_{i}.log",
        sj = "results/{sample}/rnaseq/preproc/align/bams/sj_{i}.tab"
    log:
        "logs/star_align/{sample}_inputreads_{i}.log",
    params:
        extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10"
    threads: config['threads']
    wrapper:
        "v1.26.0/bio/star/align"


rule samtools_star_merge:
    input:
        aggregate_star_align
    output:
        "results/{sample}/rnaseq/preproc/align/aln.bam",
    log:
        "logs/samtools/merge/{sample}.log",
    params:
        extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
        "v1.32.1/bio/samtools/merge"

rule samtools_star_index:
    input:
        "results/{sample}/rnaseq/preproc/align/aln.bam",
    output:
        "results/{sample}/rnaseq/preproc/align/aln.bam.bai"
    log:
        "logs/samtools/index/{sample}.log"
    params:
        extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
        "v1.31.1/bio/samtools/index"

rule samtools_postproc:
    input:
        bam="results/{sample}/rnaseq/preproc/align/aln.bam",
        idx="results/{sample}/rnaseq/preproc/align/aln.bam.bai"
    output:
        "results/{sample}/rnaseq/preproc/post/aln.bam"
    conda:
        "../envs/samtools.yml"
    log:
        "logs/samtools/postproc/{sample}.log"
    threads: 6  # more threads brings no significant increase
    shell:
        """ 
            samtools view -bh -F 4 --min-MQ {config[mapq]} {input.bam} -o - \
            | samtools sort -n -@ {threads} -m1g -O bam - -o - \
            | samtools fixmate -pcmu -O bam -@ {threads} - - \
            | samtools sort -@ {threads} -m1g -O bam - -o - \
            | samtools markdup -r -@ {threads} - {output} 
        """

