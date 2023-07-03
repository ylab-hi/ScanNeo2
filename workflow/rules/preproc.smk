if rnaseq_filetype == ".bam":
    checkpoint split_bamfile_RG:
        input:
            get_rnaseq_data
        output:
            directory("results/{sample}/rnaseq/reads/bamfiles/")
        conda:
            "../envs/samtools.yml"
        log:
            "logs/samtools/split/{sample}.logs"
        threads: 10
        shell:
            """
                mkdir -p {output}
                samtools split -@ {threads} \
                -u {output}/noRG.bam \
                -h {input} -f {output}/%!.%. {input}
            """

if rnaseq_filetype == ".bam":
    rule bam_to_fastq:
        input:
            "results/{sample}/rnaseq/reads/bamfiles/{readgroup}.bam"
        output:
            "results/{sample}/rnaseq/reads/fastqfiles/{readgroup}.fq.gz"
        conda:
            "../envs/samtools.yml"
        log:
            "logs/samtools/bam2fastq/{sample}_{readgroup}.log"
        threads: config['threads']
        shell:
            """
                samtools collate -Oun128 -@ {threads} {input} \
                | samtools fastq -OT RG -@ {threads} - | gzip -c - > {output}
            """

if rnaseq_filetype == ".bam":
    rule align_with_star:
        input:
            fq1 = "results/{sample}/rnaseq/reads/fastqfiles/{readgroup}.fq.gz",
            idx = "resources/refs/star/",
        output:
            aln = "results/{sample}/rnaseq/align/bamfiles/{readgroup}.bam",
            log = "results/{sample}/rnaseq/align/bamfiles/{readgroup}.log",
            sj = "results/{sample}/rnaseq/align/bamfiles/{readgroup}.tab"
        log:
            "logs/star_align/{sample}_{readgroup}.log"
        params:
            extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10 --outSAMattributes RG --outSAMattrRGline ID:{readgroup}"
        threads: config['threads']
        wrapper:
            "v1.26.0/bio/star/align"


if rnaseq_filetype == ".bam":
    rule merge_alignment_results:
        input:
            aggregate_alignments
        output:
            "results/{sample}/rnaseq/align/aligned.bam",
        log:
            "logs/samtools/merge/{sample}.log",
        params:
            extra="",  # optional additional parameters as string
        threads: config['threads']
        wrapper:
            "v1.32.1/bio/samtools/merge"
    

rule samtools_postproc:
    input:
        "results/{sample}/rnaseq/align/aligned.bam"
    output:
        "results/{sample}/rnaseq/align/ready.bam"
    conda:
        "../envs/samtools.yml"
    log:
        "logs/samtools/postproc/{sample}.log"
    threads: 6  # more threads brings no significant increase
    shell:
        """ 
            samtools index {input}
            samtools view -bh -F 4 --min-MQ {config[mapq]} {input} -o - \
            | samtools sort -n -@ {threads} -m1g -O bam - -o - \
            | samtools fixmate -pcmu -O bam -@ {threads} - - \
            | samtools sort -@ {threads} -m1g -O bam - -o - \
            | samtools markdup -r -@ {threads} - {output} 
        """

rule samtools_postproc_index:
    input:
        "results/{sample}/rnaseq/align/ready.bam"
    output:
        "results/{sample}/rnaseq/align/ready.bam.bai"
    log:
        "logs/samtools/index/{sample}.log"
    params:
        extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
        "v1.31.1/bio/samtools/index"


# retrieve readgroups from bam file
if rnaseq_filetype == ".bam":
    rule determine_readgroups:
        input:
            get_rnaseq_data
        output:
            "results/{sample}/rnaseq/reads/readgroups.txt"
        log:
            "logs/readgroups/{sample}.log"
        shell:
            """
                python workflow/scripts/get_readgroups.py {input} \
                {output} > {log} 2>&1
            """

rule realign:
    input:
        bam="results/{sample}/rnaseq/align/ready.bam",
        rg="results/{sample}/rnaseq/reads/readgroups.txt"
    output:
        "results/{sample}/rnaseq/align/realigned.bam"
    threads: config['threads']
    shell:
        """
            samtools collate -Oun128 {input.bam} \
            | samtools fastq -OT RG,BC - \
            | bwa mem -pt{threads} -CH <(cat {input.rg}) resources/refs/bwa/genome - \
            | samtools sort -@6 -m1g - > {output}
        """


rule realign_index:
    input:
        "results/{sample}/rnaseq/align/realigned.bam"
    output:
        "results/{sample}/rnaseq/align/realigned.bam.bai"
    log:
        "logs/samtools/index/{sample}.log"
    params:
        extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
        "v1.31.1/bio/samtools/index"



