rule trimmomatic_SE:
  input:
      unpack(get_raw_reads)
  output:
      "results/{sample}/rnaseq/preproc/reads.fq.gz"
  log:
      "logs/{sample}/trimmomatic.log"
  params:
      trimmer=["TRAILING:3"],
      extra="",
      compression_level="-9"
  threads: config['threads']
  resources:
      mem_mb=1024
  wrapper:
      "v2.1.1/bio/trimmomatic/se"

#  if rnaseq_readtype == "PE":
rule trimmomatic_PE:
  input:
    unpack(get_raw_reads)
  output:
    r1="results/{sample}/{seqtype}/reads/{replicate}_preproc_r1.fq.gz",
    r2="results/{sample}/{seqtype}/reads/{replicate}_preproc_r2.fq.gz",
    r1_unpaired="results/{sample}/{seqtype}/reads/{replicate}_preproc_r1_unpaired.fq.gz",
    r2_unpaired="results/{sample}/{seqtype}/reads/{replicate}_preproc_r2_unpaired.fq.gz"
  params:
    trimmer=[f"MINLEN:{config['preproc']['minlen']}"] 
      + [f"TRAILING:{config['preproc']['trailing']}" if config['preproc']['trailing'] is not None else ""]
      + [f"LEADING:{config['preproc']['trailing']}" if config['preproc']['leading'] is not None else ""]
      + [f"SLIDINGWINDOW:{config['preproc']['slidingwindow']['windowsize']}:{config['preproc']['slidingwindow']['quality']}" if config['preproc']['slidingwindow']['activate'] else ""]
      + [f"ILLUMINACLIP:{config['preproc']['adapters']}:2:30:10" if config['preproc']['adapters'] is not None else ""],
      extra=""
  log:
    "logs/{sample}/trimmomatic/{replicate}_{seqtype}.log"
  threads: config['threads']
  resources:
      mem_mb=1024
  wrapper:
      "v2.1.1/bio/trimmomatic/pe"

rule add_rg_fastq_PE:
  input:
#    r1="results/{sample}/{seqtype}/reads/{replicate}_preproc_r1.fq.gz",
#    r2="results/{sample}/{seqtype}/reads/{replicate}_preproc_r2.fq.gz",
    unpack(get_reads),
  output:
    r1="results/{sample}/{seqtype}/reads/{replicate}_preproc_RG_r1.fq.gz",
    r2="results/{sample}/{seqtype}/reads/{replicate}_preproc_RG_r2.fq.gz"
  message:
    "Adding read group information to fastq files"
  log:
    "logs/{sample}/add_rg/{replicate}_{seqtype}.log"
  conda:
    "../envs/basic.yml"
  shell:
    """
      bash workflow/scripts/addrgfq.sh {input.r1} > gzip -c - > {output.r1} 2> {log}      
      bash workflow/scripts/addrgfq.sh {input.r2} > gzip -c - > {output.r2} 2> {log}
    """



checkpoint splitfastq:
  input:
    unpack(get_splitfastq_input)
  output:
    directory("results/{sample}/reads/rnaseq/{replicate}/")
  log:
    "logs/{sample}/splitfastq/{replicate}.log"
  conda:
    "../envs/splitfastq.yml"
  threads: 0
  shell:
    """
      python workflow/scripts/splitfastq.py '{input}' {output} 20000000
    """

rule star_align_fastq:
  input:
    fq1 = "results/{sample}/reads/rnaseq/{replicate}/r1/reads_{i}.fq.gz",
    fq2 = "results/{sample}/reads/rnaseq/{replicate}/r2/reads_{i}.fq.gz",
    idx = "resources/refs/star/",
  output:
    aln = "results/{sample}/rnaseq/align/{replicate}/splt/reads_{i}.bam",
    log = "results/{sample}/rnaseq/align/{replicate}/splt/reads_{i}.log",
    sj = "results/{sample}/rnaseq/align/{replicate}/splt/reads_{i}.tab"
  log:
    "logs/star_align/{sample}_{replicate}_{i}.log"
  params:
      extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10 --outSAMattributes RG --outSAMattrRGline ID:noRG"
  threads: config['threads']
  wrapper:
      "v1.26.0/bio/star/align"

rule merge_alignment_results_fastq:
    input:
      aggregate_alignments_fastq
    output:
        "results/{sample}/rnaseq/align/{replicate}_aligned.bam",
    log:
        "logs/samtools/merge/{sample}_{replicate}.log",
    params:
        extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
        "v1.32.1/bio/samtools/merge"

rule star_align_pe:
    input:
        fq1 = "results/{sample}/rnaseq/reads/fastqfiles/r1/inputreads_{i}.fq.gz",
        fq2 = "results/{sample}/rnaseq/reads/fastqfiles/r2/inputreads_{i}.fq.gz",
        idx = "resources/refs/star/",
    output:
        aln = "results/{sample}/rnaseq/align/bamfiles/inputreads_{i}.bam",
        log = "results/{sample}/rnaseq/align/bamfiles/inputreads_{i}.log",
        sj = "results/{sample}/rnaseq/align/bamfiles/inputreads_{i}.tab"
    log:
        "logs/star_align/{sample}_{i}.log"
    params:
        extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10 --outSAMattributes RG --outSAMattrRGline ID:noRG"
    threads: config['threads']
    wrapper:
        "v1.26.0/bio/star/align"

rule merge_alignment_results_pe:
    input:
        aggregate_alignments_pe
    output:
        "results/{sample}/rnaseq/align/aligned.bam",
    log:
        "logs/samtools/merge/{sample}.log",
    params:
        extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
        "v1.32.1/bio/samtools/merge"


  #rule align_with_star_fq:
    #input:
        #unpack(get_align_input),
        #idx="resources/refs/star/"
    #output:
        ## see STAR manual for additional output files
        #aln="results/{sample}/rnaseq/align/aligned.bam",
        #log="logs/{sample}/star/Log1.out",
        #sj="results/{sample}/rnaseq/align/sj.out.tab"
    #log:
        #"logs/{sample}/star/Log.out"
    #params:
        #extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10 --outSAMattributes RG --outSAMattrRGline ID:xxx"
    #threads: config['threads']
    #wrapper:
        #"v2.1.1/bio/star/align"


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
        "results/{sample}/rnaseq/align/{replicate}_aligned.bam"
    output:
        "results/{sample}/rnaseq/align/{replicate}_ready.bam"
    conda:
        "../envs/samtools.yml"
    log:
        "logs/samtools/postproc/{sample}_{replicate}.log"
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
        "results/{sample}/rnaseq/align/{replicate}_ready.bam"
    output:
        "results/{sample}/rnaseq/align/{replicate}_ready.bam.bai"
    log:
        "logs/samtools/index/postproc_{sample}_{replicate}.log"
    params:
        extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
        "v1.31.1/bio/samtools/index"


# retrieve readgroups from bam file
rule get_readgroups:
    input:
      get_readgroups_input
      #"results/{sample}/rnaseq/align/{replicate}_ready.bam"
    output:
        "results/{sample}/rnaseq/reads/{replicate}_readgroups.txt"
    conda:
      "../envs/basic.yml"
    log:
        "logs/{sample}/get_readgroups/{replicate}.log"
    shell:
        """
            python workflow/scripts/get_readgroups.py '{input}' \
            {output} > {log} 2>&1
        """





rule realign:
    input:
        bam="results/{sample}/rnaseq/align/{replicate}_ready.bam",
        rg="results/{sample}/rnaseq/reads/{replicate}_readgroups.txt"
    output:
        "results/{sample}/rnaseq/align/{replicate}_realigned.bam"
    threads: config['threads']
    shell:
        """
          samtools collate -Oun128 {input.bam} \
            | samtools fastq -OT RG -@ {threads} - \
            | bwa mem -pt{threads} -CH <(cat {input.rg}) resources/refs/bwa/genome - \
            | samtools sort -@6 -m1g - > {output}
        """


rule realign_index:
    input:
        "results/{sample}/rnaseq/align/{replicate}_realigned.bam"
    output:
        "results/{sample}/rnaseq/align/{replicate}_realigned.bam.bai"
    log:
        "logs/{sample}/realign_index/{sample}_{replicate}.log"
    params:
        extra="",  # optional additional parameters as string
    threads: config['threads']
    wrapper:
        "v1.31.1/bio/samtools/index"



