if config["rnaseq-format"] == "fastq":
    rule star_align:
        input:
            fq1 = get_fqs,
            idx = "refs/star/",
        output:
            aln = "results/preproc/align/{sample}.bam",
            log = "results/preproc/align/{sample}.log",
            sj = "results/preproc/align/{sample}_SJ.out.tab",
        conda: 
            "../envs/star.yml"
        log:
            "logs/star_align_{sample}.log",
        params:
            extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --limitBAMsortRAM 11209980843"
        threads: config['threads']
        wrapper:
            "v1.26.0/bio/star/align"
    
rule filter:
    input:
        get_bams
    output:
        bam="results/preproc/post/{sample}_flt.bam",
    conda: 
        "../envs/samtools.yml"
    log:
        "logs/mapq_filter_{sample}.out",
    params:
        extra="--min-MQ {config[mapq]}",  # optional params string
        region="",  # optional region string
    threads: config['threads']
    wrapper: 
        "v1.26.0/bio/samtools/view"

rule preproc_dedup:
    input:
        bams="results/preproc/post/{sample}_flt.bam"
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="results/preproc/post/{sample}_dedup.bam",
        metrics="results/preproc/post/{sample}.metrics.txt",
    conda: 
        "../envs/picard.yml"
    log:
        "logs/dedup_bam_{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024,
    wrapper:
        "v1.26.0/bio/picard/markduplicates"

rule index:
    input:
        "results/preproc/post/{sample}_dedup.bam",
    output:
        "results/preproc/post/{sample}_dedup.bam.bai"
    log:
        "logs/samtools_index_{sample}.log"
    params:
        extra="",  # optional params string
    threads: config['threads']  # This value - 1 will be sent to -@
    wrapper:
        "v1.26.0/bio/samtools/index"
