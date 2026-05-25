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


rule sort_reads_mhcI_SE:
    input:
        "results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE.bam",
    output:
        bam="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE_sorted.bam",
    log:
        "logs/{sample}/hlatyping/sort_reads_mhcI_SE_{group}_{nartype}.log",
    conda:
        "../envs/samtools.yml"
    threads: 4
    resources:
        mem_mb=20000,
    message:
        "Sort the filtered {wildcards.nartype}seq reads for hlatyping of sample: {wildcards.sample}"
    shell:
        """
        samtools sort -n -@ {threads} -m4g {input:q} -o {output.bam:q} >{log} 2>&1
        """


checkpoint split_reads_mhcI_SE:
    input:
        fwd="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE_sorted.bam",
    output:
        directory("results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE/"),
    log:
        "logs/{sample}/hlatyping/split_reads_mhcI_SE_{group}_{nartype}.log",
    conda:
        "../envs/gatk.yml"
    threads: 1
    message:
        "Splitting filtered group: {wildcards.group} BAM files ({wildcards.nartype}seq reads) for HLA typing"
    shell:
        """
        mkdir -p results/{wildcards.sample}/hla/mhc-I/reads/{wildcards.group}_{wildcards.nartype}_flt_SE/
        gatk SplitSamByNumberOfReads \
            -I {input.fwd:q} \
            --OUTPUT {output:q} \
            --OUT_PREFIX R \
            --SPLIT_TO_N_READS 100000 \
            >{log} 2>&1
        """


rule hlatyping_mhcI_SE:
    input:
        fwd="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE/R_{no}.bam",
        rev="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_SE/R_{no}.bam",
    output:
        pdf="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE/{no}_coverage_plot.pdf",
        tsv="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE/{no}_result.tsv",
    log:
        "logs/{sample}/hlatyping/hlatyping_mhcI_SE_{group}_{nartype}_{no}.log",
    conda:
        "../envs/optitype.yml"
    # OptiType is single-threaded (ILP solver); the wrapper does not
    # parallelise across cores.
    threads: 1
    message:
        "HLA typing from splitted BAM files"
    shell:
        """
        python3 workflow/scripts/genotyping/optitype_wrapper.py \
            {wildcards.nartype} {wildcards.no} \
            results/{wildcards.sample}/hla/mhc-I/genotyping/{wildcards.group}_{wildcards.nartype}_flt_SE/ \
            {input.fwd:q} {input.rev:q} >{log} 2>&1
        """


rule combine_hlatyping_mhcI_SE:
    input:
        aggregate_mhcI_SE,
    output:
        "results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_SE.tsv",
    log:
        "logs/{sample}/hlatyping/combine_hlatyping_mhcI_SE_{group}_{nartype}.log",
    conda:
        "../envs/basic.yml"
    threads: 1
    message:
        "Combining HLA alleles from predicted optitype results from {wildcards.nartype}seq reads in group: {wildcards.group}"
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


rule sort_and_index_reads_mhcI_PE:
    input:
        "results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_{readpair}.bam",
    output:
        bam="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_{readpair}_sorted.bam",
    log:
        "logs/{sample}/hlatyping/sort_and_index_reads_mhcI_PE_{group}_{nartype}_{readpair}.log",
    conda:
        "../envs/samtools.yml"
    threads: 4
    resources:
        mem_mb=20000,
    message:
        "Sort the filtered {wildcards.nartype}seq reads for hlatyping of sample: {wildcards.sample} with readpair: {wildcards.readpair}"
    shell:
        """
        samtools sort -n -@ {threads} -m4g {input:q} -o {output.bam:q} >{log} 2>&1
        """


checkpoint split_reads_mhcI_PE:
    input:
        fwd="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_R1_sorted.bam",
        rev="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE_R2_sorted.bam",
    output:
        directory("results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE/"),
    log:
        "logs/{sample}/hlatyping/split_reads_mhcI_PE_{group}_{nartype}.log",
    conda:
        "../envs/gatk.yml"
    threads: 1
    message:
        "Splitting filtered group: {wildcards.group} BAM files ({wildcards.nartype}seq reads) for HLA typing"
    shell:
        """
        mkdir -p results/{wildcards.sample}/hla/mhc-I/reads/{wildcards.group}_{wildcards.nartype}_flt_PE/
        gatk SplitSamByNumberOfReads \
            -I {input.fwd:q} \
            --OUTPUT {output:q} \
            --OUT_PREFIX R1 \
            --SPLIT_TO_N_READS 100000 \
            >{log} 2>&1

        gatk SplitSamByNumberOfReads \
            -I {input.rev:q} \
            --OUTPUT {output:q} \
            --OUT_PREFIX R2 \
            --SPLIT_TO_N_READS 100000 \
            >>{log} 2>&1
        """


rule hlatyping_mhcI_PE:
    input:
        fwd="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE/R1_{no}.bam",
        rev="results/{sample}/hla/mhc-I/reads/{group}_{nartype}_flt_PE/R2_{no}.bam",
    output:
        pdf="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE/{no}_coverage_plot.pdf",
        tsv="results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE/{no}_result.tsv",
    log:
        "logs/{sample}/hlatyping/hlatyping_mhcI_PE_{group}_{nartype}_{no}.log",
    conda:
        "../envs/optitype.yml"
    # OptiType is single-threaded (ILP solver); the wrapper does not
    # parallelise across cores.
    threads: 1
    message:
        "HLA typing from splitted BAM files"
    shell:
        """
        python3 workflow/scripts/genotyping/optitype_wrapper.py \
            {wildcards.nartype} {wildcards.no} \
            results/{wildcards.sample}/hla/mhc-I/genotyping/{wildcards.group}_{wildcards.nartype}_flt_PE/ \
            {input.fwd:q} {input.rev:q} >{log} 2>&1
        """


rule combine_hlatyping_mhcI_PE:
    input:
        aggregate_mhcI_PE,
    output:
        "results/{sample}/hla/mhc-I/genotyping/{group}_{nartype}_flt_PE.tsv",
    log:
        "logs/{sample}/hlatyping/combine_hlatyping_mhcI_PE_{group}_{nartype}.log",
    conda:
        "../envs/basic.yml"
    threads: 1
    message:
        "Combining HLA alleles from predicted optitype results from {wildcards.nartype}seq reads in group: {wildcards.group}"
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
