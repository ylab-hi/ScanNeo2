rule fastqc_single_end:
    input:
      get_qc_input
    output:
        html="results/{sample}/{seqtype}/qualitycontrol/{group}_fastqc.html",
        zip="results/{sample}/{seqtype}/qualitycontrol/{group}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/{sample}/fastqc/{seqtype}_{group}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.2.1/bio/fastqc"

rule preproc_single_end:
  input:
    get_preproc_input,
    qc="results/{sample}/{seqtype}/qualitycontrol/{group}_fastqc.html"
  output:
    trimmed="results/{sample}/{seqtype}/reads/{group}_preproc.fq.gz",
    failed="results/{sample}/{seqtype}/reads/{group}_preproc_failed.fq.gz",
    html="results/{sample}/{seqtype}/reads/{group}_preproc_report.html",
    json="results/{sample}/{seqtype}/reads/{group}_preproc_report.json"
  log:
    "logs/{sample}/fastp/{seqtype}_{group}.log"
  params:
    adapters="",
    extra=""
  threads: config['threads']
  wrapper:
    "v2.2.1/bio/fastp"

rule fastqc_forward:
    input:
      get_qc_input_fwd
    output:
        html="results/{sample}/{seqtype}/qualitycontrol/{group}_r1_fastqc.html",
        zip="results/{sample}/{seqtype}/qualitycontrol/{group}_r1_fastqc.zip"
    params:
        extra = ""
    log:
        "logs/{sample}/fastqc/{seqtype}_{group}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.2.1/bio/fastqc"

rule fastqc_reverse:
    input:
      get_qc_input_rev
    output:
        html="results/{sample}/{seqtype}/qualitycontrol/{group}_r2_fastqc.html",
        zip="results/{sample}/{seqtype}/qualitycontrol/{group}_r2_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/{sample}/fastqc/{seqtype}_{group}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.2.1/bio/fastqc"

rule preproc_paired_end:
    input:
      unpack(get_preproc_input),
      gc1="results/{sample}/{seqtype}/qualitycontrol/{group}_r1_fastqc.html",
      gc2="results/{sample}/{seqtype}/qualitycontrol/{group}_r2_fastqc.html"
    output:
        trimmed=["results/{sample}/{seqtype}/reads/{group}_preproc_r1.fq.gz", 
                 "results/{sample}/{seqtype}/reads/{group}_preproc_r2.fq.gz"],
        unpaired1="results/{sample}/{seqtype}/reads/{group}_preproc_r1_unpaired.fq.gz",
        unpaired2="results/{sample}/{seqtype}/reads/{group}_preproc_r2_unpaired.fq.gz",
        failed="results/{sample}/{seqtype}/reads/{group}_preproc_failed.fq.gz",
        html="results/{sample}/{seqtype}/reads/{group}_preproc_report.html",
        json="results/{sample}/{seqtype}/reads/{group}_preproc_report.json",
    log:
        "logs/{sample}/fastp/{seqtype}_{group}.log"
    params:
      dapters="",
      extra=lambda wildcards: "-u 100 -e {0} -l {1} ".format(
          config['preproc']['quality'], 
          config['preproc']['minlen'])
    threads: 2
    wrapper:
        "v2.2.1/bio/fastp"


#rule trimmomatic_se:
  #input:
      #unpack(get_raw_reads)
  #output:
      #"results/{sample}/rnaseq/preproc/reads.fq.gz"
  #log:
      #"logs/{sample}/trimmomatic.log"
  #params:
      #trimmer=[f"minlen:{config['preproc']['minlen']}"] 
      #+ [f"trailing:{config['preproc']['trailing']}" if config['preproc']['trailing'] is not none else ""]
      #+ [f"leading:{config['preproc']['trailing']}" if config['preproc']['leading'] is not none else ""]
      #+ [f"slidingwindow:{config['preproc']['slidingwindow']['windowsize']}:{config['preproc']['slidingwindow']['quality']}" if config['preproc']['slidingwindow']['activate'] else ""]
      #+ [f"illuminaclip:{config['preproc']['adapters']}:2:30:10" if config['preproc']['adapters'] is not none else ""],
      #extra="",
      #compression_level="-9"
  #threads: config['threads']
  #resources:
      #mem_mb=1024
  #wrapper:
      #"v2.1.1/bio/trimmomatic/se"

##  if rnaseq_readtype == "pe":
#rule trimmomatic_pe:
  #input:
    #unpack(get_raw_reads)
  #output:
    #r1="results/{sample}/{seqtype}/reads/{replicate}_preproc_r1.fq.gz",
    #r2="results/{sample}/{seqtype}/reads/{replicate}_preproc_r2.fq.gz",
    #r1_unpaired="results/{sample}/{seqtype}/reads/{replicate}_preproc_r1_unpaired.fq.gz",
    #r2_unpaired="results/{sample}/{seqtype}/reads/{replicate}_preproc_r2_unpaired.fq.gz"
  #params:
    #trimmer=[f"minlen:{config['preproc']['minlen']}"] 
      #+ [f"trailing:{config['preproc']['trailing']}" if config['preproc']['trailing'] is not none else ""]
      #+ [f"leading:{config['preproc']['trailing']}" if config['preproc']['leading'] is not none else ""]
      #+ [f"slidingwindow:{config['preproc']['slidingwindow']['windowsize']}:{config['preproc']['slidingwindow']['quality']}" if config['preproc']['slidingwindow']['activate'] else ""]
      #+ [f"illuminaclip:{config['preproc']['adapters']}:2:30:10" if config['preproc']['adapters'] is not none else ""],
      #extra=""
  #log:
    #"logs/{sample}/trimmomatic/{replicate}_{seqtype}.log"
  #threads: config['threads']
  #resources:
      #mem_mb=1024
  #wrapper:
      #"v2.1.1/bio/trimmomatic/pe"

#rule add_rg_fastq_pe:
  #input:
    #unpack(get_reads),
  #output:
    #r1="results/{sample}/{seqtype}/reads/{replicate}_preproc_rg_r1.fq.gz",
    #r2="results/{sample}/{seqtype}/reads/{replicate}_preproc_rg_r2.fq.gz"
  #message:
    #"adding read group information to fastq files"
  #log:
    #"logs/{sample}/add_rg/{replicate}_{seqtype}.log"
  #conda:
    #"../envs/basic.yml"
  #shell:
    #"""
      #bash workflow/scripts/addrgfq.sh {input.r1} > gzip -c - > {output.r1} 2> {log}      
      #bash workflow/scripts/addrgfq.sh {input.r2} > gzip -c - > {output.r2} 2> {log}
    #"""


#checkpoint splitfastq:
  #input:
    #unpack(get_splitfastq_input)
  #output:
    #directory("results/{sample}/reads/rnaseq/{replicate}/")
  #log:
    #"logs/{sample}/splitfastq/{replicate}.log"
  #conda:
    #"../envs/splitfastq.yml"
  #threads: 0
  #shell:
    #"""
      #python workflow/scripts/splitfastq.py '{input}' {output} 20000000
    #"""

#rule star_align_fastq:
  #input:
    #fq1 = "results/{sample}/reads/rnaseq/{replicate}/r1/reads_{i}.fq.gz",
    #fq2 = "results/{sample}/reads/rnaseq/{replicate}/r2/reads_{i}.fq.gz",
    #idx = "resources/refs/star/",
  #output:
    #aln = "results/{sample}/rnaseq/align/{replicate}/splt/reads_{i}.bam",
    #log = "results/{sample}/rnaseq/align/{replicate}/splt/reads_{i}.log",
    #sj = "results/{sample}/rnaseq/align/{replicate}/splt/reads_{i}.tab"
  #log:
    #"logs/star_align/{sample}_{replicate}_{i}.log"
  #params:
      #extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10 --outSAMattributes RG --outSAMattrRGline ID:noRG"
  #threads: config['threads']
  #wrapper:
      #"v1.26.0/bio/star/align"

#rule merge_alignment_results_fastq:
    #input:
      #aggregate_alignments_fastq
    #output:
        #"results/{sample}/rnaseq/align/{replicate}_aligned.bam",
    #log:
        #"logs/samtools/merge/{sample}_{replicate}.log",
    #params:
        #extra="",  # optional additional parameters as string
    #threads: config['threads']
    #wrapper:
        #"v1.32.1/bio/samtools/merge"

#rule star_align_pe:
    #input:
        #fq1 = "results/{sample}/rnaseq/reads/fastqfiles/r1/inputreads_{i}.fq.gz",
        #fq2 = "results/{sample}/rnaseq/reads/fastqfiles/r2/inputreads_{i}.fq.gz",
        #idx = "resources/refs/star/",
    #output:
        #aln = "results/{sample}/rnaseq/align/bamfiles/inputreads_{i}.bam",
        #log = "results/{sample}/rnaseq/align/bamfiles/inputreads_{i}.log",
        #sj = "results/{sample}/rnaseq/align/bamfiles/inputreads_{i}.tab"
    #log:
        #"logs/star_align/{sample}_{i}.log"
    #params:
        #extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10 --outSAMattributes RG --outSAMattrRGline ID:noRG"
    #threads: config['threads']
    #wrapper:
        #"v1.26.0/bio/star/align"

#rule merge_alignment_results_pe:
    #input:
        #aggregate_alignments_pe
    #output:
        #"results/{sample}/rnaseq/align/aligned.bam",
    #log:
        #"logs/samtools/merge/{sample}.log",
    #params:
        #extra="",  # optional additional parameters as string
    #threads: config['threads']
    #wrapper:
        #"v1.32.1/bio/samtools/merge"


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


#if rnaseq_filetype == ".bam":
    #checkpoint split_bamfile_RG:
        #input:
            #get_rnaseq_data
        #output:
            #directory("results/{sample}/rnaseq/reads/bamfiles/")
        #conda:
            #"../envs/samtools.yml"
        #log:
            #"logs/samtools/split/{sample}.logs"
        #threads: 10
        #shell:
            #"""
                #mkdir -p {output}
                #samtools split -@ {threads} \
                #-u {output}/noRG.bam \
                #-h {input} -f {output}/%!.%. {input}
            #"""

    #rule bam_to_fastq:
        #input:
            #"results/{sample}/rnaseq/reads/bamfiles/{readgroup}.bam"
        #output:
            #"results/{sample}/rnaseq/reads/fastqfiles/{readgroup}.fq.gz"
        #conda:
            #"../envs/samtools.yml"
        #log:
            #"logs/samtools/bam2fastq/{sample}_{readgroup}.log"
        #threads: config['threads']
        #shell:
            #"""
                #samtools collate -Oun128 -@ {threads} {input} \
                #| samtools fastq -OT RG -@ {threads} - | gzip -c - > {output}
            #"""

    #rule align_with_star:
        #input:
            #fq1 = "results/{sample}/rnaseq/reads/fastqfiles/{readgroup}.fq.gz",
            #idx = "resources/refs/star/",
        #output:
            #aln = "results/{sample}/rnaseq/align/bamfiles/{readgroup}.bam",
            #log = "results/{sample}/rnaseq/align/bamfiles/{readgroup}.log",
            #sj = "results/{sample}/rnaseq/align/bamfiles/{readgroup}.tab"
        #log:
            #"logs/star_align/{sample}_{readgroup}.log"
        #params:
            #extra="--outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimOutType WithinBAM HardClip --genomeSAindexNbases 10 --outSAMattributes RG --outSAMattrRGline ID:{readgroup}"
        #threads: config['threads']
        #wrapper:
            #"v1.26.0/bio/star/align"

#if rnaseq_filetype == ".bam":
    #rule merge_alignment_results:
        #input:
            #aggregate_alignments
        #output:
            #"results/{sample}/rnaseq/align/aligned.bam",
        #log:
            #"logs/samtools/merge/{sample}.log",
        #params:
            #extra="",  # optional additional parameters as string
        #threads: config['threads']
        #wrapper:
            #"v1.32.1/bio/samtools/merge"
    


#rule samtools_postproc_index:
    #input:
        #"results/{sample}/rnaseq/align/{replicate}_ready.bam"
    #output:
        #"results/{sample}/rnaseq/align/{replicate}_ready.bam.bai"
    #log:
        #"logs/samtools/index/postproc_{sample}_{replicate}.log"
    #params:
        #extra="",  # optional additional parameters as string
    #threads: config['threads']
    #wrapper:
        #"v1.31.1/bio/samtools/index"





