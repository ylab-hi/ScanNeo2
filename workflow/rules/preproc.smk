rule fastqc_single_end:
    input:
      get_qc_input
    output:
        html="results/{sample}/{seqtype}/qualitycontrol/{group}_raw.html",
        zip="results/{sample}/{seqtype}/qualitycontrol/{group}_raw.zip"
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
    sample=get_preproc_input,
    qc="results/{sample}/{seqtype}/qualitycontrol/{group}_raw.html"
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
        html="results/{sample}/{seqtype}/qualitycontrol/{group}_R1_fastqc_raw.html",
        zip="results/{sample}/{seqtype}/qualitycontrol/{group}_R1_fastqc_raw.zip"
    params:
        extra = ""
    log:
        "logs/{sample}/fastqc/{seqtype}_{group}_fwd_raw.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.2.1/bio/fastqc"

rule fastqc_reverse:
    input:
      get_qc_input_rev
    output:
        html="results/{sample}/{seqtype}/qualitycontrol/{group}_R2_fastqc_raw.html",
        zip="results/{sample}/{seqtype}/qualitycontrol/{group}_R2_fastqc_raw.zip"
    params:
        extra = "--quiet"
    log:
        "logs/{sample}/fastqc/{seqtype}_{group}_fwd_raw.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.2.1/bio/fastqc"

rule preproc_paired_end:
    input:
      unpack(get_preproc_input),
      qc1="results/{sample}/{seqtype}/qualitycontrol/{group}_R1_fastqc_raw.html",
      qc2="results/{sample}/{seqtype}/qualitycontrol/{group}_R2_fastqc_raw.html"
    output:
        trimmed=["results/{sample}/{seqtype}/reads/{group}_R1_preproc.fq.gz", 
                 "results/{sample}/{seqtype}/reads/{group}_R2_preproc.fq.gz"],
        unpaired1="results/{sample}/{seqtype}/reads/{group}_R1_preproc_unpaired.fq.gz",
        unpaired2="results/{sample}/{seqtype}/reads/{group}_R2_preproc_unpaired.fq.gz",
        failed="results/{sample}/{seqtype}/reads/{group}_preproc_failed.fq.gz",
        html="results/{sample}/{seqtype}/reads/{group}_preproc_report.html",
        json="results/{sample}/{seqtype}/reads/{group}_preproc_report.json",
    log:
        "logs/{sample}/fastp/{seqtype}_{group}.log"
    params:
      dapters="",
      extra=lambda wildcards: "-u 100 -e {0} -l {1} {2} ".format(
          config['preproc']['qual'], 
          config['preproc']['minlen'],
          "-3 --cut_tail_window_size {} cut_tail_mean_quality {}".format(
              config['preproc']['slidingwindow']['wsize'],
              config['preproc']['slidingwindow']['wqual'] if config['preproc']['slidingwindow']['activate'] else ""
              )
      )
    threads: config['threads']
    wrapper:
        "v2.2.1/bio/fastp"

# fastqc after pre-processing (forward)
rule fastqc_forward_after:
    input:
      "results/{sample}/{seqtype}/reads/{group}_R1_preproc.fq.gz"
    output:
        html="results/{sample}/{seqtype}/qualitycontrol/{group}_R1_fastqc.html",
        zip="results/{sample}/{seqtype}/qualitycontrol/{group}_R1_fastqc.zip"
    params:
        extra = ""
    log:
        "logs/{sample}/fastqc/{seqtype}_{group}_fwd.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.2.1/bio/fastqc"

# fastqc after pre-processing (reverse)
rule fastqc_reverse_after:
    input:
      "results/{sample}/{seqtype}/reads/{group}_R2_preproc.fq.gz"
    output:
        html="results/{sample}/{seqtype}/qualitycontrol/{group}_R2_fastqc.html",
        zip="results/{sample}/{seqtype}/qualitycontrol/{group}_R2_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/{sample}/fastqc/{seqtype}_{group}_rev.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v2.2.1/bio/fastqc"


