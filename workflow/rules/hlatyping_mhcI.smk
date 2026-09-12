####### CLASS I HLA GENOTYPING ###########


######### single-end reads #########
rule filter_reads_mhcI_SE:
    input:
        reads=get_input_filtering_hlatyping_SE,
        panel=multiext(
            "resources/hla/yara_index/{nartype}",
            ".lf.drp",
            ".lf.drs",
            ".lf.drv",
            ".lf.pst",
            ".rid.concat",
            ".rid.limits",
            ".sa.ind",
            ".sa.len",
            ".sa.val",
            ".txt.concat",
            ".txt.limits",
            ".txt.size",
        ),
    output:
        reads="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE.bam",
    log:
        "logs/{sample}/hlatyping/filter_reads_mhcI_SE_{group}_{nartype}.log",
    conda:
        "../envs/yara.yml"
    threads: config["threads"]
    message:
        "Filter the {wildcards.nartype} reads of group:{wildcards.group} of sample:{wildcards.sample} against the HLA panel"
    shell:
        """
        (
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara_index/{wildcards.nartype} \
                {input.reads:q} \
                | samtools view -h -F 4 -b1 - -o {output.reads:q}
        ) >{log} 2>&1
        """


# HLA typing runs OptiType once over the whole HLA-filtered read set. OptiType
# solves for a single best 6-allele genotype from all reads at once, so the
# reads must NOT be split: combine_optitype_results unions the alleles across
# inputs, and splitting would let per-piece sampling noise inflate the call
# set. The wrapper coordinate-sorts and indexes the BAM before running.
rule hlatyping_mhcI_SE:
    input:
        fwd="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE.bam",
        rev="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE.bam",
    output:
        pdf="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE_coverage_plot.pdf",
        tsv="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE_result.tsv",
    log:
        "logs/{sample}/hlatyping/hlatyping_mhcI_SE_{group}_{nartype}.log",
    conda:
        "../envs/optitype.yml"
    # OptiType is single-threaded (ILP solver); its razers3 hit matrix over the
    # whole read set is the memory driver.
    threads: 1
    resources:
        mem_mb=64000,
        runtime=360,
    message:
        "HLA typing (OptiType) of {wildcards.nartype}seq reads in group: {wildcards.group}"
    shell:
        """
        python3 workflow/scripts/genotyping/optitype_wrapper.py \
            {wildcards.nartype} {wildcards.group}_{wildcards.nartype}_flt_SE \
            results/{wildcards.sample}/hla/mhc-I/genotyping/ \
            {input.fwd:q} {input.rev:q} >{log} 2>&1
        """


rule combine_hlatyping_mhcI_SE:
    input:
        "results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE_result.tsv",
    output:
        "results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE.tsv",
    log:
        "logs/{sample}/hlatyping/combine_hlatyping_mhcI_SE_{group}_{nartype}.log",
    conda:
        "../envs/basic.yml"
    threads: 1
    message:
        "Reformatting OptiType alleles from {wildcards.nartype}seq reads in group: {wildcards.group}"
    shell:
        """
        python3 workflow/scripts/genotyping/combine_optitype_results.py \
            '{input}' {wildcards.group} {output} >{log} 2>&1
        """


############# paired-end reads ###########
rule filter_reads_mhcI_PE:
    input:
        reads=get_input_filtering_hlatyping_PE,
        panel=multiext(
            "resources/hla/yara_index/{nartype}",
            ".lf.drp",
            ".lf.drs",
            ".lf.drv",
            ".lf.pst",
            ".rid.concat",
            ".rid.limits",
            ".sa.ind",
            ".sa.len",
            ".sa.val",
            ".txt.concat",
            ".txt.limits",
            ".txt.size",
        ),
    output:
        reads="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_{readpair}.bam",
    log:
        "logs/{sample}/hlatyping/filter_reads_mhcI_PE_{group}_{nartype}_{readpair}.log",
    conda:
        "../envs/yara.yml"
    threads: config["threads"]
    message:
        "Filter the {wildcards.nartype} reads of group:{wildcards.group} of sample:{wildcards.sample} against the HLA panel"
    shell:
        """
        (
            yara_mapper -t {threads} -e 3 -f bam -u resources/hla/yara_index/{wildcards.nartype} \
                {input.reads:q} \
                | samtools view -h -F 4 -b1 - -o {output.reads:q}
        ) >{log} 2>&1
        """


# Single OptiType call over the whole R1/R2 filtered set (see hlatyping_mhcI_SE).
rule hlatyping_mhcI_PE:
    input:
        fwd="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_R1.bam",
        rev="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_R2.bam",
    output:
        pdf="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE_coverage_plot.pdf",
        tsv="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE_result.tsv",
    log:
        "logs/{sample}/hlatyping/hlatyping_mhcI_PE_{group}_{nartype}.log",
    conda:
        "../envs/optitype.yml"
    threads: 1
    resources:
        mem_mb=64000,
        runtime=360,
    message:
        "HLA typing (OptiType) of {wildcards.nartype}seq reads in group: {wildcards.group}"
    shell:
        """
        python3 workflow/scripts/genotyping/optitype_wrapper.py \
            {wildcards.nartype} {wildcards.group}_{wildcards.nartype}_flt_PE \
            results/{wildcards.sample}/hla/mhc-I/genotyping/ \
            {input.fwd:q} {input.rev:q} >{log} 2>&1
        """


rule combine_hlatyping_mhcI_PE:
    input:
        "results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE_result.tsv",
    output:
        "results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE.tsv",
    log:
        "logs/{sample}/hlatyping/combine_hlatyping_mhcI_PE_{group}_{nartype}.log",
    conda:
        "../envs/basic.yml"
    threads: 1
    message:
        "Reformatting OptiType alleles from {wildcards.nartype}seq reads in group: {wildcards.group}"
    shell:
        """
        python3 workflow/scripts/genotyping/combine_optitype_results.py \
            '{input}' {wildcards.group} {output} >{log} 2>&1
        """


rule combine_all_mhcI_alleles:
    input:
        get_all_mhcI_alleles,
    output:
        "results/{sample}/hla/mhc-I.tsv",
    log:
        "logs/{sample}/hlatyping/combine_all_mhcI_alleles.log",
    conda:
        "../envs/basic.yml"
    threads: 1
    message:
        "Combining HLA alleles from different sources (e.g., predicted and user-defined alleles)"
    shell:
        """
        python workflow/scripts/genotyping/combine_all_alleles.py \
            '{input}' mhc-I {output} >{log} 2>&1
        """
