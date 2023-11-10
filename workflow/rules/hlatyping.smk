####### CLASS I HLA GENOTYPING ###########
rule get_hla_panel:
  output:
    dna="resources/hla/hla_ref_dna.fasta",
    rna="resources/hla/hla_ref_rna.fasta"
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
        dna="resources/hla/hla_ref_dna.fasta",
        rna="resources/hla/hla_ref_rna.fasta"
    output:
        dna=multiext("resources/hla/yara/idx/dna", 
            ".lf.drp", ".lf.drs", ".lf.drv", 
            ".lf.pst", ".rid.concat", ".rid.limits",
            ".sa.ind", ".sa.len", ".sa.val", 
            ".txt.concat", ".txt.limits", ".txt.size"),
        rna=multiext("resources/hla/yara/idx/rna", 
            ".lf.drp", ".lf.drs", ".lf.drv", 
            ".lf.pst", ".rid.concat", ".rid.limits",
            ".sa.ind", ".sa.len", ".sa.val", 
            ".txt.concat", ".txt.limits", ".txt.size")
    log:
        "logs/yara_indexer.log"
    conda:
        "../envs/yara.yml"
    shell:
        """
        yara_indexer -o resources/hla/yara/idx/dna {input.dna} > {log}
        yara_indexer -o resources/hla/yara/idx/rna {input.rna} >> {log}
        """

# map reads against reference DNA
if config['data']['dnaseq'] is not None:
  if config['hlatyping']['mode'] == 'BOTH' or config['hlatyping']['mode'] == 'DNA':
    if config['data']['dnaseq_readtype'] == 'SE':
      rule get_hla_filtering_input_single_DNA:
        input:
          get_hla_flt_dna_se,
          dna=multiext("resources/hla/yara/idx/dna", 
            ".lf.drp", ".lf.drs", ".lf.drv", 
            ".lf.pst", ".rid.concat", ".rid.limits",
            ".sa.ind", ".sa.len", ".sa.val", 
            ".txt.concat", ".txt.limits", ".txt.size"),
        output:
          "results/{sample}/hla/{group}_flt_DNA.bam"
        message:
          "Mapping HLA reads against reference"
        log:
          "logs/{sample}/{group}_hla_reads_filtering_DNA"
        conda:
          "../envs/yara.yml"
        shell:
          """
            samtools fastq {input[0]} > results/{wildcards.sample}/hla/{wildcards.group}_DNA.fastq 
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/dna \
              results/{wildcards.sample}/hla/{wildcards.group}_DNA.fastq \
              | samtools view -h -F 4 -b1 - | samtools sort - -o {output}
            samtools index {output}
          """

if config['data']['dnaseq'] is not None:
  if config['hlatyping']['mode'] == 'BOTH' or config['hlatyping']['mode'] == 'DNA':
    if config['data']['dnaseq_readtype'] == 'PE':
      rule get_hla_filtering_input_paired_DNA:
        input:
          unpack(get_hla_flt_dna_pe),
          dna=multiext("resources/hla/yara/idx/dna", 
            ".lf.drp", ".lf.drs", ".lf.drv", 
            ".lf.pst", ".rid.concat", ".rid.limits",
            ".sa.ind", ".sa.len", ".sa.val", 
            ".txt.concat", ".txt.limits", ".txt.size"),
        output:
          fwd="results/{sample}/hla/{group}_R1_flt_DNA.bam",
          rev="results/{sample}/hla/{group}_R2_flt_DNA.bam",
        message:
          "Mapping HLA reads against reference"
        log:
          "logs/{sample}/{group}_hla_reads_filtering_dna"
        conda:
          "../envs/yara.yml"
        threads: config['threads']
        shell:
          """
            samtools fastq {input.r1} > results/{wildcards.sample}/hla/{wildcards.group}_flt_R1_DNA.fastq
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/dna \
                results/{wildcards.sample}/hla/{wildcards.group}_flt_R1_DNA.fastq \
                | samtools view -h -F 4 -b1 - | samtools sort - -o {output.fwd}
            samtools index {output.fwd}
            samtools fastq {input.r2} > results/{wildcards.sample}/hla/{wildcards.group}_flt_R2_DNA.fastq
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/dna \
                results/{wildcards.sample}/hla/{wildcards.group}_flt_R2_DNA.fastq \
                | samtools view -h -F 4 -b1 - | samtools sort - -o {output.rev}
            samtools index {output.rev}
          """

# map reads against reference RNA
if config['data']['rnaseq'] is not None:
  if config['hlatyping']['mode'] == 'BOTH' or config['hlatyping']['mode'] == 'RNA':
    if config['data']['rnaseq_readtype'] == 'SE':
      rule get_hla_filtering_input_single_RNA:
        input:
          get_hla_flt_rna_se,
          dna=multiext("resources/hla/yara/idx/rna", 
            ".lf.drp", ".lf.drs", ".lf.drv", 
            ".lf.pst", ".rid.concat", ".rid.limits",
            ".sa.ind", ".sa.len", ".sa.val", 
            ".txt.concat", ".txt.limits", ".txt.size"),
        output:
          "results/{sample}/hla/{group}_flt_RNA.bam"
        message:
          "Mapping HLA reads against reference"
        log:
          "logs/{sample}/{group}_hla_reads_filtering_rna.log"
        threads: config['threads']
        conda:
          "../envs/yara.yml"
        shell:
          """
            samtools fastq {input[0]} > results/{wildcards.sample}/hla/{wildcards.group}.fastq 
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/rna \
                results/{wildcards.sample}/hla/{wildcards.group}.fastq \
                | samtools view -h -F 4 -b1 - | samtools sort - -o {output}
            samtools index {output}
          """

if config['data']['rnaseq'] is not None:
  if config['hlatyping']['mode'] == 'BOTH' or config['hlatyping']['mode'] == 'RNA':
    if config['data']['rnaseq_readtype'] == 'PE':
      rule get_hla_filtering_input_paired_RNA:
        input:
          unpack(get_hla_flt_rna_pe),
          dna=multiext("resources/hla/yara/idx/rna", 
            ".lf.drp", ".lf.drs", ".lf.drv", 
            ".lf.pst", ".rid.concat", ".rid.limits",
            ".sa.ind", ".sa.len", ".sa.val", 
            ".txt.concat", ".txt.limits", ".txt.size"),
        output:
          fwd="results/{sample}/hla/{group}_R1_flt_RNA.bam",
          rev="results/{sample}/hla/{group}_R2_flt_RNA.bam",
        message:
          "Mapping HLA reads against reference"
        log:
          "logs/{sample}/{group}_hla_reads_filtering_RNA"
        conda:
          "../envs/yara.yml"
        threads: config['threads']
        shell:
          """
            samtools fastq {input.r1} > results/{wildcards.sample}/hla/{wildcards.group}_r1.fastq
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/rna \
                results/{wildcards.sample}/hla/{wildcards.group}_r1.fastq \
                | samtools view -h -F 4 -b1 - | samtools sort - -o {output.fwd}
            samtools index -m1g {output.fwd}
            samtools fastq {input.r2} > results/{wildcards.sample}/hla/{wildcards.group}_r2.fastq
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/rna \
                results/{wildcards.sample}/hla/{wildcards.group}_r2.fastq \
                | samtools view -h -F 4 -b1 - | samtools sort - -o {output.rev}
            samtools index -m1g {output.rev}
          """

checkpoint split_bam_single:
  input:
    "results/{sample}/hla/{group}_flt_{type}.bam",
  output:
    directory("results/{sample}/hla/{group}_flt_{type}_splitted"),
  message:
    "Splitting filtered BAM files for HLA typing"
  log:
    "logs/{sample}/hla/split_bam_{group}_{type}.log"
  conda:
    "../envs/gatk.yml"
  threads: 1
  shell:
    """
      mkdir -p results/{wildcards.sample}/hla/{wildcards.group}_flt_{wildcards.type}_splitted/
      gatk SplitSamByNumberOfReads \
          -I {input[0]} \
          --OUTPUT {output} \
          --OUT_PREFIX R \
          --SPLIT_TO_N_READS 100000
    """


rule hla_genotyping_single:
  input:
    "results/{sample}/hla/{group}_flt_{type}_splitted/R_{no}.bam",
  output:
    pdf="results/{sample}/hla/{group}_flt_{type}_typed/{no}_coverage_plot.pdf",
    tsv="results/{sample}/hla/{group}_flt_{type}_typed/{no}_result.tsv"
  message:
    "HLA typing from splitted BAM files"
  log:
    "logs/{sample}/optitype/{group}_{type}_{no}_call.log"
  conda:
    "../envs/optitype.yml"
  threads: config['threads']
  shell:
    """
      samtools index {input}
      if [ "{wildcards.type}" == "RNA" ]; then
        OptiTypePipeline.py --input {input} \
            --outdir results/{wildcards.sample}/hla/{wildcards.group}_flt_{wildcards.type}_typed/ \
            --prefix {wildcards.no} --rna -v > {log}
      elif [ "{wildcards.type}" == "DNA" ]; then
        OptiTypePipeline.py --input {input} \
            --outdir results/{wildcards.sample}/hla/{wildcards.group}_flt_{wildcards.type}_typed/ \
            --prefix {wildcards.no} --dna -v > {log}
      fi
    """

rule combine_genotyping_single:
  input:
    aggregate_genotyping_single
  output:
    "results/{sample}/hla/alleles/classI_{group}_{type}_SE.tsv"
  log:
    "logs/{sample}/optitype/{group}_{type}_call.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      python3 workflow/scripts/merge_classI_alleles.py \
          '{input}' {output}
    """ 

checkpoint split_bam_paired:
  input:
    fwd="results/{sample}/hla/{group}_R1_flt_{type}.bam",
    rev="results/{sample}/hla/{group}_R2_flt_{type}.bam"
  output:
    directory("results/{sample}/hla/{group}_flt_{type}_splitted"),
  message:
    "Splitting filtered BAM files for HLA typing"
  log:
    "logs/{sample}/hla/split_bam_{group}_{type}.log"
  conda:
    "../envs/gatk.yml"
  threads: 1
  shell:
    """
      mkdir -p results/{wildcards.sample}/hla/{wildcards.group}_flt_{wildcards.type}_splitted/
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


rule hla_genotyping_paired:
  input:
    fwd="results/{sample}/hla/{group}_flt_{type}_splitted/R1_{no}.bam",
    rev="results/{sample}/hla/{group}_flt_{type}_splitted/R2_{no}.bam"
  output:
    pdf="results/{sample}/hla/{group}_flt_{type}_typed/{no}_coverage_plot.pdf",
    tsv="results/{sample}/hla/{group}_flt_{type}_typed/{no}_result.tsv"
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
      if [ "{wildcards.type}" == "RNA" ]; then
        OptiTypePipeline.py --input {input.fwd} {input.rev} \
            --outdir results/{wildcards.sample}/hla/{wildcards.group}_flt_{wildcards.type}_typed/ \
            --prefix {wildcards.no} --rna -v > {log}
      elif [ "{wildcards.type}" == "DNA" ]; then
        OptiTypePipeline.py --input {input.fwd} {input.rev} \
            --outdir results/{wildcards.sample}/hla/{wildcards.group}_flt_{wildcards.type}_typed/ \
            --prefix {wildcards.no} --dna -v > {log}
      fi
    """

rule combine_genotyping_paired:
  input:
    aggregate_genotyping_paired
  output:
    "results/{sample}/hla/alleles/classI_{group}_{type}_PE.tsv"
  log:
    "logs/{sample}/optitype/{group}_{type}_call.log"
  shell:
    """
      python3 workflow/scripts/merge_classI_alleles.py \
          '{input}' {output}
    """ 


rule merge_allels:
  input:
    get_alleles
  output:
    "results/{sample}/hla/classI_alleles.tsv",
  message:
    "Merging HLA alleles from different sources"
  log:
    "logs/{sample}/optitype/merge_classI_alleles.log"
  conda:
    "../envs/basic.yml"
  threads: 1
  shell:
    """
      cat {input} | sort | uniq > {output}
    """

####### CLASS II HLA GENOTYPING ###########
#rule classII_get_reads_DNA_SE:
  #input:
    #get_hla_flt_dna_se
  #output:
    #"results/{sample}/hla/classII/{group}_DNA.fastq"
  #message:
    #"Extracting DNA reads for HLA classII typing"
  #log:
    #"logs/{sample}/optitype/{group}_DNA_call.log"
  #conda:
    #"../envs/basic.yml"
  #threads: 1
  #shell:
    #"""
      
    #"""


#def get_hla_flt_dna_se(wildcards):

