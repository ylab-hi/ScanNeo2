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
            directory("results/{sample}/rnaseq/reads/fastqfiles/{readgroup}.fq.gz")
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
        "results/{sample}/rnaseq/preproc/aligned.bam"
    output:
        "results/{sample}/rnaseq/preproc/ready.bam"
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
        "results/{sample}/rnaseq/preproc/ready.bam"
    output:
        "results/{sample}/rnaseq/preproc/ready.bam.bai"
    log:
        "logs/samtools/index/{sample}.log"
    params:
        extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
        "v1.31.1/bio/samtools/index"


# create 
#rule get_readgroups:
#    input:
#        get_rnaseq_data
#    output:
#        star="results/{sample}/rnaseq/reads/readgroups_star.txt"
#        bwa="results/{sample}/rnaseq/reads/readgroups_bwa.txt"
#    conda:
#        "../envs/readgroups.yml"






# if the case of specified bam files, convert to fastq
#if rnaseq_filetype == ".bam":
#    rule bamfile_to_fastq:
#        input:
#            get_rnaseq_data
#        output:
#            "results/{sample}/rnaseq/reads/inputreads.fq.gz"
#        conda:
#            "../envs/samtools.yml"
#        log:
#            "logs/samtools/bam2fastq/{sample}.log"
#        threads: config['threads']
#        shell:
#            """
#                samtools collate -Oun128 -@ {threads} {input} \
#                | samtools fastq -OT RG -@ {threads} - | gzip -c - > {output}
#            """


# create interleaved format for fastq file
# if rnaseq_filetype == ".fastq":
# rule interleaved_fastq:
#   input:
#       get_rnaseq_data
#   output:
#       "results/{sample}/rnaseq/reads/inputreads.fq"
#   conda:
#       "../envs/samtools.yml"
#   log:
#       "logs/samtools/bam2fastq/{sample}.log"
#   shell:
#       """
#       """


# create
#checkpoint split_fastqfile:
#    input: 
#        "results/{sample}/rnaseq/reads/inputreads.fq.gz"
#    output:
#        directory("results/{sample}/rnaseq/reads/fastqfiles/")
#    conda:
#        "../envs/splitfastq.yml"
#    log:
#        "logs/splitfastq/{sample}.log"
#    shell:
#        """
#            python workflow/scripts/splitfastq.py {input} \
#            results/{wildcards.sample}/rnaseq/reads/fastqfiles/ \
#            10000000 > {log} 2>&1
#        """


# further split the bamfiles by readgroup
#rule star_align:
#    input:
#        fq1 = "results/{sample}/rnaseq/reads/fastqfiles/inputreads_{file}.fq.gz",
#        idx = "resources/refs/star/",
#    output:
#        "results/{sample}/rnaseq/align/bamfiles/inputreads_{file}.bam",
#    conda:
#        "../envs/star.yml"
#    log:
#        "logs/star_align/{sample}_{file}.log"
#    threads:
#        config['threads']
#    shell:
#        """
#            python workflow/scripts/star_align.py {input.fq1} \
#            {threads} {input.idx} {output} > {log} 2>&1
#        """
        


#if rnaseq_filetype == ".bam":
#    checkpoint split_bamfile:
#        input:
#            get_rnaseq_data
#        output:
#            directory("results/{sample}/rnaseq/reads/bamfiles/")
#        conda:
#            "../envs/picard.yml"
#        log:
#            "logs/picard/split/{sample}.logs"
#        shell:
#            """
#                mkdir -p {output}
#                picard SplitSamByNumberOfReads --INPUT {input} \
#                --OUTPUT {output} --OUT_PREFIX inputreads \
#                --SPLIT_TO_N_READS 10000000 > {log} 2>&1
#           """










#        input:
#            get_bams,
#        output:
#            directory("results/{sample}/rnaseq/reads/bamfiles/")
#        conda:
#            "../envs/samtools.yml"
#        log:
#            "logs/samtools/split/{sample}.logs"
#        threads: 10
#        shell:
#            """
#                mkdir -p {output}
#                samtools split -@ {threads} \
#                -u {output}/noRG.bam \
#                -h {input} -f {output}/%!.%. {input}
#            """


#rule star_align:
#
#    input:
#        fq1="results/{sample}/rnaseq/reads/fqs/inputreads_{i}.fq.gz",
#        idx="resources/refs/star/",






        








    










#    checkpoint bamfile_split:
#        input:
#            get_bams,
#        output:
#            directory("results/{sample}/rnaseq/reads/bamfiles/")
#        conda:
#            "../envs/samtools.yml"
#        log:
#            "logs/samtools/split/{sample}.logs"
#        threads: 10
#        shell:
#            """
#                mkdir -p {output}
#                samtools split -@ {threads} \
#                -u {output}/noRG.bam \
#                -h {input} -f {output}/%!.%. {input}
#            """


#if rnaseq_filetype == ".bam":
#    rule bamfile_readgroups:
#        input:
#            get_bams,
#        output:
#            "results/{sample}/rnaseq/reads/readgroups.txt"
#        conda:
#            "../envs/samtools.yml"
#        shell:
#            "samtools view -H {input} | grep ^@RG > {output}"



#rule bamfile_define_readgroups:
#i#    input:
  #      "results/{sample}/rnaseq/reads/readgroups.txt"
  #  output:
  #      "results/{sample}/rnaseq/reads/readgroups_modified.txt"
  #i  run:
   #     shell(touch {output})
        #print("read groups file")
#        with open("{input}","r") as f:
#            for line in f:
#                print(line)



# reads read group from input file or dummy
#if rnaseq_filetype == ".bam":
#    rule define_readgroups:
#        input:
#            get_bams, 
#        output: 
#            "results/{sample}/rnaseq/reads/readgroups.txt",
#            "results/{sample}/rnaseq/reads/readgroups_modified.txt"
#        conda:
#            "../envs/samtools.yml"
#        run:
#            """
#                shell(samtools view -H {input} | grep ^@RG > {output[0]})
#                with open("{output[0]}","r") as f:
#                    for line in f:
#                        print(line)
#                shell(touch {output[1]})
#
#            """





# split bam files into read groups
# align the reads with STAR


























        
#if rnaseq_filetype == ".bam":
#    rule bamfile_split_and_align:
#        input:
#            "results/{sample}/rnaseq/reads/{file}.bam"
#        output:
#            "results/{sample}/rnaseq/align/bamfiles/{file}.bam"
#        conda:
#            "../envs/align.yml"
#        threads: config['threads']
#        shell:
#            """
#                samtools collate -Oun128 {input} \
#                | samtools fastq -OT RG,BC - \
#                | bwa mem -pt{threads} -CH <(samtools view -H {input}|grep ^@RG) resources/refs/bwa/genome - \
#                | samtools sort -@6 -m1g - > {output}
#            """


# create rules for .fastq files (when input is fastq)
# if rnaseq_filetype == ".fastq" or rnaseq_filetype == ".fq":




#checkpoint split_aligned_bamfile:
#    input:
#        "results/{sample}/rnaseq/postproc/aln.bam"
#    output:
#        directory("results/{sample}/rnaseq/postproc/splitted/")
#    log:
#        "logs/splitbam/{sample}.log"
#    shell:
#        """
#            python workflow/scripts/splitbam.py {input} 2000000 {output} 
#        """

# align with STAR
#rule star_align:
#    input:
#        fq1 = "results/{sample}/rnaseq/postproc/splitted/aln_{i}.bam",
#        idx = "resources/refs/star/",
#    output:
#        aln = "results/{sample}/rnaseq/postproc/splitted/realn_{i}.bam",
#        log = "results/{sample}/rnaseq/postproc/splitted/realn_{i}.log",
#        sj = "results/{sample}/rnaseq/preproc/realn_{i}.tab" 
#    log:
#        "logs/star_align/{sample}_{i}.log"
#    params:
#        extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10"
#    threads: config['threads']
#    wrapper:
#        "v1.26.0/bio/star/align"

#rule realigned_bamfiles_merge: 
#    input:
#        aggregate_star_alignments
#    output:
#        "results/{sample}/rnaseq/postproc/realn.bam",
#    log:
#        "logs/samtools/merge/{sample}.log",
#    params:
#        extra="",  # optional additional parameters as string
#    threads: config['threads']
#    wrapper:
#        "v1.32.1/bio/samtools/merge"


#rule samtools_realign_index:
#    input:
#        "results/{sample}/rnaseq/postproc/realn.bam"
#    output:
#        "results/{sample}/rnaseq/postproc/realn.bam.bai"
#    log:
#        "logs/samtools/index/{sample}.log"
#    params:
#        extra="",  # optional additional parameters as string
#    threads: config['threads']
#    wrapper:
#        "v1.31.1/bio/samtools/index"


# if the case of specified bam files, convert to fastq
#rule bam_to_fastq:
#    input:
#        "results/{sample}/rnaseq/preproc/bams/{file}.bam"
#    output:
#        "results/{sample}/rnaseq/preproc/pre/fqs/{file}.fq.gz"
#    conda:
#        "../envs/samtools.yml"
#    log:
#        "logs/samtools/bam2fastq/{sample}.log"
#    threads: config['threads']
#    shell:
#        """
#            samtools collate -Oun128 -@ {threads} {input} \
#            | samtools fastq -OT RG,BC -@ {threads} - | gzip -c - > {output}
#        """


#rule align_bwa:
#    input:
#        "results/{sample}/rnaseq/preproc/pre/fqs/{file}.fq.gz",
#    output:
#        "results{sample}/rnaseq/align/bams/{file}.bam"
#    conda:
#        "../envs/realign.yml"
#    log:
#        "logs/bwa/{sample}.log"
#    threads: config['threads']
#    shell:
#        """
#            bwa mem -pt{threads} -C
#        """
    











#if rnaseq_filetype == ".bam":
#    rule bam_to_fastq:
#        input:
#            get_bams,
#        output:
#            "results/{sample}/rnaseq/preproc/pre/inputreads.fq.gz"
#        conda:
#            "../envs/samtools.yml"
#        log:
#            "logs/samtools/bam2fastq/{sample}.log"
#        threads: config['threads']
#        shell:
#            """
#                samtools collate -Oun128 -@ {threads} {input} \
#                | samtools fastq -OT RG,BC -@ {threads} - | gzip -c - > {output}
#            """

#checkpoint split_fastq:
#    input:
#        "results/{sample}/rnaseq/preproc/pre/inputreads.fq.gz"
#    output:
#        directory("results/{sample}/rnaseq/preproc/pre/fqs/")
#    log:
#        "logs/split_fastq/{sample}.log"
#    conda:
#        "../envs/splitfastq.yml"
#    shell:
#        """
#            python workflow/scripts/splitfastq.py {input} \
#            results/{wildcards.sample}/rnaseq/preproc/pre/ \
#            10000000 > {log} 2>&1
#        """

#extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --limitBAMsortRAM 11209980843"
#rule star_align:
#    input:
#        fq1 = "results/{sample}/rnaseq/preproc/pre/fqs/inputreads_{i}.fq.gz",
#        idx = "resources/refs/star/",
#    output:
#        aln = "results/{sample}/rnaseq/preproc/align/bams/aln_{i}.bam",
#        log = "results/{sample}/rnaseq/preproc/align/bams/star.align_{i}.log",
#        sj = "results/{sample}/rnaseq/preproc/align/bams/sj_{i}.tab"
#    log:
#        "logs/star_align/{sample}_inputreads_{i}.log",
#    params:
#        extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10"
#    threads: config['threads']
#    wrapper:
#        "v1.26.0/bio/star/align"

#rule samtools_star_merge:
#    input:
#        aggregate_star_align
#    output:
#        "results/{sample}/rnaseq/preproc/align/aln.bam",
#    log:
#        "logs/samtools/merge/{sample}.log",
#    params:
#        extra="",  # optional additional parameters as string
#    threads: config['threads']
#    wrapper:
#        "v1.32.1/bio/samtools/merge"

#rule samtools_star_index:
#    input:
#        "results/{sample}/rnaseq/preproc/align/aln.bam",
#    output:
#        "results/{sample}/rnaseq/preproc/align/aln.bam.bai"
#    log:
#        "logs/samtools/index/{sample}.log"
#    params:
#        extra="",  # optional additional parameters as string
#    threads: config['threads']
#    wrapper:
#        "v1.31.1/bio/samtools/index"


