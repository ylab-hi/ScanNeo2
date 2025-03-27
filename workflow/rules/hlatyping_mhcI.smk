####### CLASS I HLA GENOTYPING ###########

######### single-end reads #########
rule filter_reads_mhcI_SE:
  input:
    reads=get_input_filtering_hlatyping_SE,
    panel=multiext("resources/hla/yara_index/{nartype}",
        ".lf.drp", ".lf.drs", ".lf.drv",
        ".lf.pst", ".rid.concat", ".rid.limits",
        ".sa.ind", ".sa.len", ".sa.val", 
        ".txt.concat", ".txt.limits", ".txt.size"),
  output:
    reads="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE.bam",
  message:
    "Filter the {wildcards.nartype} reads of group:{wildcards.group} of sample:{wildcards.sample} against the HLA panel"
  log:
    "logs/{sample}/hlatyping/filter_reads_mhcI_SE_{group}_{nartype}.log"
  conda:
    "../envs/yara.yml"
  threads: config['threads']
  shell:
    """
      yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara_index/{wildcards.nartype} \
          {input.reads} | samtools view -h -F 4 -b1 - -o {output.reads} > {log} 2>&1
    """

rule sort_and_index_reads_mhcI_SE:
  input:
    "results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE.bam",
  output:
    bam="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE_sorted.bam",
    idx="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE_sorted.bam.bai"
  message:
    "Sort the filtered {wildcards.nartype}seq reads for hlatyping of sample: {wildcards.sample}"
  threads: 4
  resources:
    mem_mb=20000
  log:
    "logs/{sample}/hlatyping/sort_and_index_reads_mhcI_PE_{group}_{nartype}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      samtools sort -@ {threads} -m4g {input} -o {output.bam} > {log} 2>&1
      samtools index {output.bam} >> {log} 2>&1
    """

checkpoint split_reads_mhcI_SE:
  input:
    fwd="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE_sorted.bam",
    fwd_idx="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE_sorted.bam.bai",
  output:
    directory("results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE/")
  message:
    "Splitting filtered group: {wildcards.group} BAM files ({wildcards.nartype}seq reads) for HLA typing"
  log:
    "logs/{sample}/hlatyping/split_bam_mhcI_SE_{group}_{nartype}.log"
  conda:
    "../envs/gatk.yml"
  threads: 1
  shell:
    """
      mkdir -p results/{wildcards.sample}/hla/mhc-I/reads/{wildcards.group}_{wildcards.nartype}_flt_SE/
      gatk SplitSamByNumberOfReads \
          -I {input.fwd} \
          --OUTPUT {output} \
          --OUT_PREFIX R \
          --SPLIT_TO_N_READS 100000
    """

rule hlatyping_mhcI_SE:  
  input:
    fwd="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE/R_{no}.bam",
    rev="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE/R_{no}.bam",
  output:
    pdf="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE/{no}_coverage_plot.pdf",
    tsv="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE/{no}_result.tsv"
  message:
    "HLA typing from splitted BAM files"
  log:
    "logs/{sample}/hlatyping/hlatyping_mhcI_SE_{group}_{nartype}_{no}.log"
  conda:
    "../envs/optitype.yml"
  threads: 64
  shell:
    """
      python3 workflow/scripts/genotyping/optitype_wrapper.py \
          '{input.fwd} {input.rev}' {wildcards.nartype} {wildcards.no} \
          results/{wildcards.sample}/hla/mhc-I/genotyping/{wildcards.group}_{wildcards.nartype}_flt_SE/ > {log} 2>&1
    """

rule combine_hlatyping_mhcI_SE:
  input:
    aggregate_mhcI_SE,
  output:
    "results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE.tsv",
  message:
    "Combining HLA alleles from predicted optitype results from {wildcards.nartype}seq reads in group: {wildcards.group}"
  log:
    "logs/{sample}/hlatyping/combine_optitype_mhcI_PE_{group}_{nartype}.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python3 workflow/scripts/genotyping/combine_optitype_results.py \
          '{input}' {wildcards.group} {output}
    """

############# paired-end reads ###########
rule filter_reads_mhcI_PE:
  input:
    reads=get_input_filtering_hlatyping_PE,
    panel=multiext("resources/hla/yara_index/{nartype}",
        ".lf.drp", ".lf.drs", ".lf.drv",
        ".lf.pst", ".rid.concat", ".rid.limits",
        ".sa.ind", ".sa.len", ".sa.val", 
        ".txt.concat", ".txt.limits", ".txt.size"),
  output:
    reads="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_{readpair}.bam",
  message:
    "Filter the {wildcards.nartype} reads of group:{wildcards.group} of sample:{wildcards.sample} against the HLA panel"
  log:
    "logs/{sample}/hlatyping/filter_reads_mhcI_PE_{group}_{nartype}_{readpair}.log"
  conda:
    "../envs/yara.yml"
  threads: config['threads']
  shell:
    """
      yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara_index/{wildcards.nartype} \
          {input.reads} | samtools view -h -F 4 -b1 - -o {output.reads} > {log} 2>&1
    """

rule sort_and_index_reads_mhcI_PE:
  input:
    "results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_{readpair}.bam",
  output:
    bam="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_{readpair}_sorted.bam",
    idx="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_{readpair}_sorted.bam.bai"
  message:
    "Sort the filtered {wildcards.nartype}seq reads for hlatyping of sample: {wildcards.sample} with readpair: {wildcards.readpair}"
  threads: 4
  resources:
    mem_mb=20000
  log:
    "logs/{sample}/hlatyping/sort_and_index_reads_mhcI_PE_{group}_{nartype}_{readpair}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      samtools sort -@ {threads} -m4g {input} -o {output.bam} > {log} 2>&1
      samtools index {output.bam} >> {log} 2>&1
    """

checkpoint split_reads_mhcI_PE:
  input:
    fwd="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_R1_sorted.bam",
    fwd_idx="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_R1_sorted.bam.bai",
    rev="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_R2_sorted.bam",
    rev_idx="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_R2_sorted.bam.bai"
  output:
    directory("results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE/")
  message:
    "Splitting filtered group: {wildcards.group} BAM files ({wildcards.nartype}seq reads) for HLA typing"
  log:
    "logs/{sample}/hlatyping/split_bam_mhcI_PE_{group}_{nartype}.log"
  conda:
    "../envs/gatk.yml"
  threads: 1
  shell:
    """
      mkdir -p results/{wildcards.sample}/hla/mhc-I/reads/{wildcards.group}_{wildcards.nartype}_flt_PE/
      gatk SplitSamByNumberOfReads \
          -I {input.fwd} \
          --OUTPUT {output} \
          --OUT_PREFIX R1 \
          --SPLIT_TO_N_READS 100000
      
      gatk SplitSamByNumberOfReads \
          -I {input.rev} \
          --OUTPUT {output} \
          --OUT_PREFIX R2 \
          --SPLIT_TO_N_READS 100000
    """

rule hlatyping_mhcI_PE:  
  input:
    fwd="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE/R1_{no}.bam",
    rev="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE/R2_{no}.bam",
  output:
    pdf="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE/{no}_coverage_plot.pdf",
    tsv="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE/{no}_result.tsv"
  message:
    "HLA typing from splitted BAM files"
  log:
    "logs/{sample}/hlatyping/hlatyping_mhcI_PE_{group}_{nartype}_{no}.log"
  conda:
    "../envs/optitype.yml"
  threads: 64
  shell:
    """
      python3 workflow/scripts/genotyping/optitype_wrapper.py \
          '{input.fwd} {input.rev}' {wildcards.nartype} {wildcards.no} \
          results/{wildcards.sample}/hla/mhc-I/genotyping/{wildcards.group}_{wildcards.nartype}_flt_PE/ > {log} 2>&1
    """

rule combine_hlatyping_mhcI_PE:
  input:
    aggregate_mhcI_PE,
  output:
    "results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE.tsv",
  message:
    "Combining HLA alleles from predicted optitype results from {wildcards.nartype}seq reads in group: {wildcards.group}"
  log:
    "logs/{sample}/hlatyping/combine_optitype_mhcI_PE_{group}_{nartype}.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python3 workflow/scripts/genotyping/combine_optitype_results.py \
          '{input}' {wildcards.group} {output}
    """

rule combine_all_mhcI_alleles:
  input:
    get_all_mhcI_alleles
  output:
    "results/{sample}/hla/mhc-I.tsv"
  message:
    "Combining HLA alleles from different sources (e.g., predicted and user-defined alleles)"
  log:
    "logs/{sample}/genotyping/combine_all_mhc-I.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python workflow/scripts/genotyping/combine_all_alleles.py \
          '{input}' mhc-I {output}
    """

