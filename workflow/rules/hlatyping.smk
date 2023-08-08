rule get_hla_panel:
    output:
        dna="resources/hla/hla_ref_dna.fasta",
        rna="resources/hla/hla_ref_rna.fasta"
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
          "results/{sample}/hla/{group}_flt_dna.bam"
        message:
          "Mapping HLA reads against reference"
        log:
          "logs/{sample}/{group}_hla_reads_filtering_dna"
        conda:
          "../envs/yara.yml"
        shell:
          """
            yara_mapper -e 3 -f bam -u resources/hla/yara/idx/dna {input[0]} \
                | samtools view -h -F 4 -b1 - > {output}
          """

if config['data']['dnaseq'] is not None:
  if config['hlatyping']['mode'] == 'BOTH' or config['hlatyping']['mode'] == 'DNA':
    if config['data']['dnaseq_readtype'] == 'PE':
      rule get_hla_filtering_input_paired_DNA:
        input:
          get_hla_flt_dna_pe,
          dna=multiext("resources/hla/yara/idx/dna", 
            ".lf.drp", ".lf.drs", ".lf.drv", 
            ".lf.pst", ".rid.concat", ".rid.limits",
            ".sa.ind", ".sa.len", ".sa.val", 
            ".txt.concat", ".txt.limits", ".txt.size"),
        output:
          "results/{sample}/hla/{group}_flt_r1_dna.bam",
          "results/{sample}/hla/{group}_flt_r2_dna.bam",
        message:
          "Mapping HLA reads against reference"
        log:
          "logs/{sample}/{group}_hla_reads_filtering_dna"
        conda:
          "../envs/yara.yml"
        threads: config['threads']
        shell:
          """
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/dna {input[0]} \
                | samtools view -h -F 4 -b1 - | samtools sort - -o {output[0]} > {log}
            samtools index {output[0]}
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/dna {input[1]} \
                | samtools view -h -F 4 -b1 - | samtools sort - -o {output[1]} >> {log}
            samtools index {output[1]}
          """


# map reads against reference RNA
if config['data']['rnaseq'] is not None:
  if config['hlatyping']['mode'] == 'BOTH' or config['hlatyping']['mode'] == 'RNA':
    if config['data']['rnaseq_readtype'] == 'SE':
      rule get_hla_filtering_input_single_RNA:
        input:
          get_hla_flt_rna_se,
          dna=multiext("resources/hla/yara/idx/dna", 
            ".lf.drp", ".lf.drs", ".lf.drv", 
            ".lf.pst", ".rid.concat", ".rid.limits",
            ".sa.ind", ".sa.len", ".sa.val", 
            ".txt.concat", ".txt.limits", ".txt.size"),
        output:
          "results/{sample}/hla/{group}_flt_rna.bam"
        message:
          "Mapping HLA reads against reference"
        log:
          "logs/{sample}/{group}_hla_reads_filtering_rna"
        threads: config['threads']
        conda:
          "../envs/yara.yml"
        shell:
          """
            samtools fastq {input[0]} > results/{wildcards.sample}/hla/{wildcards.group}_rna.fastq 
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/rna \
                results/{wildcards.sample}/hla/{wildcards.group}_rna.fastq \
                | samtools view -h -F 4 -b1 - | samtools sort - -o {output}
            samtools index {output}
          """

if config['data']['rnaseq'] is not None:
  if config['hlatyping']['mode'] == 'BOTH' or config['hlatyping']['mode'] == 'RNA':
    if config['data']['dnaseq_readtype'] == 'PE':
      rule get_hla_filtering_input_paired_RNA:
        input:
          get_hla_flt_rna_pe,
          dna=multiext("resources/hla/yara/idx/dna", 
            ".lf.drp", ".lf.drs", ".lf.drv", 
            ".lf.pst", ".rid.concat", ".rid.limits",
            ".sa.ind", ".sa.len", ".sa.val", 
            ".txt.concat", ".txt.limits", ".txt.size"),
        output:
          "results/{sample}/hla/{group}_flt_r1_rna.bam",
          "results/{sample}/hla/{group}_flt_r2_rna.bam",
        message:
          "Mapping HLA reads against reference"
        log:
          "logs/{sample}/{group}_hla_reads_filtering_rna"
        conda:
          "../envs/yara.yml"
        threads: config['threads']
        shell:
          """
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/rna {input[0]} \
                | samtools view -h -F 4 -b1 - > {output[0]}
            samtools index {output[0]}
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara/idx/rna {input[1]} \
                | samtools view -h -F 4 -b1 - > {output[1]}
            samtools index {output[1]}
          """

rule hla_genotyping_DNA:
  input:
    unpack(hlatyping_input_DNA),
  output:
    pdf="results/{sample}/hla/{group}_dna_coverage_plot.pdf",
    tsv="results/{sample}/hla/{group}_dna_result.tsv",
  log:
    "logs/{sample}/optitype/{group}_call.log"
  conda:
    "../envs/optitype.yml"
  shell:
    """
      OptiTypePipeline.py --input {input} \
          --outdir results/{wildcards.sample}/hla/ \
          --prefix {wildcards.group}_dna --dna > {log}
    """

rule hla_genotyping_RNA:
  input:
    unpack(hlatyping_input_RNA),
  output:
    pdf="results/{sample}/hla/{group}_rna_coverage_plot.pdf",
    tsv="results/{sample}/hla/{group}_rna_result.tsv",
  log:
    "logs/{sample}/optitype/{group}_call.log"
  conda:
    "../envs/optitype.yml"
  shell:
    """
      OptiTypePipeline.py --input {input} \
          --outdir results/{wildcards.sample}/hla/ \
          --prefix {wildcards.group}_dna --dna > {log}
    """


### merge






