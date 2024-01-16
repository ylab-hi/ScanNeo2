####### CLASS I HLA GENOTYPING ###########
rule get_hla_panel:
  output:
    dna="resources/hla/hla_ref_DNA.fa",
    rna="resources/hla/hla_ref_RNA.fa"
  message:
    "Downloading HLA reference panels"
  conda:
    "../envs/basic.yml"
  log:
    "logs/hla_panel.log"
  shell:
    """
      curl -o {output.dna} https://raw.githubusercontent.com/FRED-2/OptiType/v1.3.5/data/hla_reference_dna.fasta
      curl -o {output.rna} https://raw.githubusercontent.com/FRED-2/OptiType/v1.3.5/data/hla_reference_rna.fasta
    """

rule index_hla_panel:
    input:
        fa="resources/hla/hla_ref_{type}.fa",
    output:
        panel=multiext("resources/hla/yara_index/{type}", 
            ".lf.drp", ".lf.drs", ".lf.drv", 
            ".lf.pst", ".rid.concat", ".rid.limits",
            ".sa.ind", ".sa.len", ".sa.val", 
            ".txt.concat", ".txt.limits", ".txt.size"),
    log:
        "logs/yara_indexer_{type}.log"
    conda:
        "../envs/yara.yml"
    shell:
      """
          yara_indexer -o resources/hla/yara_index/{wildcards.type} {input.fa} > {log}
      """

######### get input reads ########
rule get_reads_hlatyping_SE:
  input:
    get_input_hlatyping_SE
  output:
    "results/{sample}/hla/reads/{group}_{type}_SE.fq"
  message:
    "Retrieve reads for HLA genotyping by converting the alignments ({wildcards.type}) of group:{wildcards.group} from sample:{wildcards.sample} to reads in FASTQ format"
  log:
    "logs/{sample}/hla/get_hla_input_reads_{group}_{type}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      samtools fastq {input} > {output}
    """

rule get_reads_hlatyping_PE:
  input:
    unpack(get_input_hlatyping_PE)
  output:
    fwd="results/{sample}/hla/reads/{group}_{type}_PE_{readtype}.fq",
  message:
    "Retrieve paired-end reads ({wildcards.type}) for HLA genotyping - sample:{wildcards.sample} group:{wildcards.group}"
  log:
    "logs/{sample}/hla/get_reads_hlatyping_PE_{group}_{type}_{readtype}.log"
  conda:
    "../envs/samtools.yml"
  shell:
    """
      samtools fastq {input.fwd} > {output.fwd}
    """

######### single-end reads #########
rule filter_reads_mhcI_SE:
  input:
    reads="results/{sample}/hla/reads/{group}_{type}_SE.fq",
    panel=multiext("resources/hla/yara_index/{type}", 
        ".lf.drp", ".lf.drs", ".lf.drv", 
        ".lf.pst", ".rid.concat", ".rid.limits",
        ".sa.ind", ".sa.len", ".sa.val", 
        ".txt.concat", ".txt.limits", ".txt.size"),
  output:
    "results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_SE.bam"
  message:
    "Filter the {wildcards.type} reads of group:{wildcards.group} of sample:{wildcards.sample} against the HLA panel"
  log:
    "logs/{sample}/{group}_hla_reads_filtering_{type}"
  threads:
    config["threads"]
  conda:
    "../envs/yara.yml"
  shell:
    """
      if [ "{wildcards.type}" == "DNA" ]; then
        yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara_index/{wildcards.type} \
          {input.reads} | samtools view -h -F 4 -b1 - | samtools sort - -o {output} > {log}
      elif [ "{wildcards.type}" == "RNA" ]; then
        yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara_index/{wildcards.type} \
            {input.reads} | samtools view -h -F 4 -b1 - | samtools sort - -o {output} > {log}
      fi
    """

rule index_reads_mhcI_SE:
  input:
    "results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_SE.bam"
  output:
    "results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_SE.bam.bai"
  log:
    "logs/{sample}/{group}_hla_reads_filtering_{type}_index"
  params:
      extra="",  # optional params string
  threads: 4  # This value - 1 will be sent to -@
  wrapper:
    "v2.3.0/bio/samtools/index"

checkpoint split_reads_mhcI_SE:
  input:
    reads="results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_SE.bam",
    index="results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_SE.bam.bai"
  output:
    directory("results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_SE/")
  message:
    "Splitting filtered BAM files for HLA typing"
  log:
    "logs/{sample}/hla/split_bam_{group}_{type}.log"
  conda:
    "../envs/gatk.yml"
  threads: 1
  shell:
    """
      mkdir -p results/{wildcards.sample}/hla/mhc-I/reads/{wildcards.group}_{wildcards.type}_flt_SE/
      gatk SplitSamByNumberOfReads \
          -I {input.reads} \
          --OUTPUT {output} \
          --OUT_PREFIX R \
          --SPLIT_TO_N_READS 100000
    """

rule hlatyping_mhcI_SE:  
  input:
    reads="results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_SE/R_{no}.bam",
  output:
    pdf="results/{sample}/hla/mhc-I/genotyping/{group}_{type}_flt_SE/{no}_coverage_plot.pdf",
    tsv="results/{sample}/hla/mhc-I/genotyping/{group}_{type}_flt_SE/{no}_result.tsv"
  message:
    "HLA typing from splitted BAM files"
  log:
    "logs/{sample}/optitype/{group}_{type}_{no}_call.log"
  conda:
    "../envs/optitype.yml"
  threads: config['threads']
  shell:
    """
      samtools index {input.reads}
      if [ "{wildcards.type}" == "DNA" ]; then
        OptiTypePipeline.py --input {input.reads} \
            --outdir results/{wildcards.sample}/hla/mhc-I/genotyping/{wildcards.group}_{wildcards.type}_flt_SE/ \
            --prefix {wildcards.no} --dna -v > {log}
      elif [ "{wildcards.type}" == "RNA" ]; then
        OptiTypePipeline.py --input {input.reads} \
            --outdir results/{wildcards.sample}/hla/mhc-I/genotyping/{wildcards.group}_{wildcards.type}_flt_SE/ \
            --prefix {wildcards.no} --rna -v > {log}
      fi
    """

rule combine_mhcI_SE:
  input:
    aggregate_mhcI_SE
  output:
    "results/{sample}/hla/mhc-I/genotyping/{group}_{type}_SE.tsv"
  log:
    "logs/{sample}/optitype/{group}_{type}_call.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python3 workflow/scripts/combine_optitype_results.py \
          '{input}' {output}
    """

############# paired-end reads ###########
rule filter_reads_mhcI_PE:
  input:
    reads="results/{sample}/hla/reads/{group}_{type}_PE_{readtype}.fq",
    panel=multiext("resources/hla/yara_index/{type}",
        ".lf.drp", ".lf.drs", ".lf.drv",
        ".lf.pst", ".rid.concat", ".rid.limits",
        ".sa.ind", ".sa.len", ".sa.val", 
        ".txt.concat", ".txt.limits", ".txt.size"),
  output:
    reads="results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_PE_{readtype}.bam",
  message:
    "Filter the {wildcards.type} reads of group:{wildcards.group} of sample:{wildcards.sample} against the HLA panel"
  log:
    "logs/{sample}/{group}_hla_reads_filtering_{type}_{readtype}.log"
  conda:
    "../envs/yara.yml"
  threads: config['threads']
  shell:
    """
      yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara_index/{wildcards.type} \
          {input.reads} | samtools view -h -F 4 -b1 - | samtools sort -O bam - -o {output.reads} > {log}
    """

rule index_reads_mhcI_PE:
  input:
    "results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_PE_{readtype}.bam"
  output:
    "results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_PE_{readtype}.bam.bai"
  log:
    "logs/{sample}/{group}_hla_reads_filtering_{type}_index_{readtype}"
  params:
      extra="",  # optional params string
  threads: 4  # This value - 1 will be sent to -@
  wrapper:
    "v2.3.0/bio/samtools/index"

checkpoint split_reads_mhcI_PE:
  input:
    fwd="results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_PE_R1.bam",
    fwd_idx="results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_PE_R1.bam.bai",
    rev="results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_PE_R2.bam",
    rev_idx="results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_PE_R2.bam.bai"
  output:
    directory("results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_PE/")
  message:
    "Splitting filtered BAM files for HLA typing"
  log:
    "logs/{sample}/hla/split_bam_{group}_{type}.log"
  conda:
    "../envs/gatk.yml"
  threads: 1
  shell:
    """
      mkdir -p results/{wildcards.sample}/hla/mhc-I/reads/{wildcards.group}_{wildcards.type}_flt_PE/
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
    fwd="results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_PE/R1_{no}.bam",
    rev="results/{sample}/hla/mhc-I/reads/{group}_{type}_flt_PE/R2_{no}.bam",
  output:
    pdf="results/{sample}/hla/mhc-I/genotyping/{group}_{type}_flt_PE/{no}_coverage_plot.pdf",
    tsv="results/{sample}/hla/mhc-I/genotyping/{group}_{type}_flt_PE/{no}_result.tsv"
  message:
    "HLA typing from splitted BAM files"
  log:
    "logs/{sample}/optitype/{group}_{type}_{no}_call.log"
  conda:
    "../envs/optitype.yml"
  threads: config['threads']
  shell:
    """
      samtools index {input.fwd}
      samtools index {input.rev}
      if [ "{wildcards.type}" == "DNA" ]; then
        OptiTypePipeline.py --input {input.fwd} {input.rev} \
            --outdir results/{wildcards.sample}/hla/mhc-I/genotyping/{wildcards.group}_{wildcards.type}_flt_PE/ \
            --prefix {wildcards.no} --dna -v > {log}
      elif [ "{wildcards.type}" == "RNA" ]; then
        OptiTypePipeline.py --input {input.fwd} {input.rev} \
            --outdir results/{wildcards.sample}/hla/mhc-I/genotyping/{wildcards.group}_{wildcards.type}_flt_PE/ \
            --prefix {wildcards.no} --rna -v > {log}
      fi
    """

rule combine_mhcI_PE:
  input:
    aggregate_mhcI_PE
  output:
    "results/{sample}/hla/mhc-I/genotyping/{group}_{type}_PE.tsv"
  log:
    "logs/{sample}/optitype/{group}_{type}_call.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python3 workflow/scripts/combine_optitype_results.py \
          '{input}' {output}
    """

rule merge_predicted_mhcI_allels:
  input:
    get_predicted_mhcI_alleles
  output:
    "results/{sample}/hla/genotyping/mhc-I.tsv",
  message:
    "Merging HLA alleles from different sources"
  log:
    "logs/{sample}/optitype/merge_predicted_mhc-I.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python workflow/scripts/genotyping/merge_predicted_mhcI.py \
          '{input}' {output}
    """

rule combine_all_mhcI_alleles:
  input:
    get_all_mhcI_alleles:
  output:
    "results/{sample}/hla/mhc-I.tsv"
  message:
    "Combining HLA alleles from different sources"
  log:
    "logs/{sample}/genotyping/combine_all_mhc-I.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python workflow/scripts/genotyping/combine_all_alleles.py \
          '{input}' {output} > {log} 2>&1\
    """
    
######### MHC-II HLA GENOTYPING ###########
rule filter_reads_mhcII_PE:
  input:
    sample=["results/{sample}/hla/reads/{group}_{type}_PE_R1.fq", "results/{sample}/hla/reads/{group}_{type}_PE_R2.fq"],
    idx=multiext(
      "resources/hla/bowtie2_index",
      ".1.bt2",
      ".2.bt2",
      ".3.bt2",
      ".4.bt2",
      ".rev.1.bt2",
      ".rev.2.bt2",
    ),
  output:
    "results/{sample}/hla/mhc-II/reads/{group}_{type}_hla.bam"
  log:
    "logs/bowtie2/{sample}/{group}_{type}.log",
  params:
    extra="",  # optional parameters
  threads: config['threads']  # Use at least two threads
  wrapper:
    "v2.11.1/bio/bowtie2/align"

rule bam2fastq_reads_mhcII:
  input:
    "results/{sample}/hla/mhc-II/reads/{group}_{type}_hla.bam"
  output:
    fwd="results/{sample}/hla/mhc-II/reads/{group}_{type}_hla_R1.fastq",
    rev="results/{sample}/hla/mhc-II/reads/{group}_{type}_hla_R2.fastq",
  log:
    "logs/{sample}/{group}_{type}_hla_sam.log"
  conda:
    "../envs/samtools.yml"
  threads: 1
  shell:
    """
      samtools fastq -F 4 {input} -1 {output.fwd} -2 {output.rev} -s /dev/null
    """

rule hlatyping_mhcII:
  input:
    fwd="results/{sample}/hla/mhc-II/reads/{group}_{type}_hla_R1.fastq",
    rev="results/{sample}/hla/mhc-II/reads/{group}_{type}_hla_R2.fastq"
  output:
    "results/{sample}/hla/mhc-II/genotyping/{group}_{type}/result/{group}_{type}_final.result.txt"
  log:
    "logs/{sample}/hla/{group}_{type}_hlahd.log"
  conda:
    "../envs/hlahd.yml"
  params:
    freqdata=f"""-f {config['hlatyping']['freqdata']}""",
    split=f"""{config['hlatyping']['split']}""",
    dic=f"""{config['hlatyping']['dict']}"""
  threads: config['threads']
  shell:
    """
      hlahd.sh \
          -t {threads} \
          -m 100 \
          -c 0.95 \
          {params.freqdata} \
          {input.fwd} {input.rev} \
          {params.split} {params.dic} \
          {wildcards.group}_{wildcards.type} \
          results/{wildcards.sample}/hla/mhc-II/genotyping/
    """

#rule merge_mhcII_allels:
  #input:
    #get_mhcII_alleles
  #output:
    #"results/{sample}/hla/mhc-I.tsv",
  #message:
    #"Merging HLA alleles from different sources"
  #log:
    #"logs/{sample}/optitype/merge_classI_alleles.log"
  #conda:
    #"../envs/basic.yml"
  #threads: 1
  #shell:
    #"""
      #python workflow/scripts/merge_mhcI_alleles.py \
          #'{input}' {output}
    #"""



    











