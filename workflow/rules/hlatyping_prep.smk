############################ hlatyping_prep #####################################
rule get_mhcI_hla_panel:
    output:
        dna="resources/hla/hla_ref_DNA.fa",
        rna="resources/hla/hla_ref_RNA.fa",
    log:
        "logs/ref/get_mhcI_hla_panel.log",
    conda:
        "../envs/basic.yml"
    message:
        "Downloading MHC-I HLA reference panels"
    shell:
        """
        (
            curl --fail -o {output.dna} https://raw.githubusercontent.com/FRED-2/OptiType/v1.3.5/data/hla_reference_dna.fasta
            curl --fail -o {output.rna} https://raw.githubusercontent.com/FRED-2/OptiType/v1.3.5/data/hla_reference_rna.fasta
        ) >{log} 2>&1
        """


rule index_mhcI_hla_panel:
    input:
        fa="resources/hla/hla_ref_{nartype}.fa",
    output:
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
    log:
        "logs/ref/index_mhcI_hla_panel_{nartype}.log",
    conda:
        "../envs/yara.yml"
    message:
        "Index MHC-I HLA reference panel ({wildcards.nartype})"
    shell:
        """
        mkdir -p resources/hla/yara_index/
        yara_indexer \
            -o resources/hla/yara_index/{wildcards.nartype} \
            {input.fa} >{log} 2>&1
        """


rule get_reads_hlatyping_BAM:
    input:
        reads=get_input_reads_hlatyping_BAM,
    output:
        "results/{sample}/hla/reads/{group}_{nartype}_BAM.fq",
    log:
        "logs/{sample}/hlatyping/get_reads_hlatyping_BAM_{group}_{nartype}.log",
    conda:
        "../envs/samtools.yml"
    threads: 4
    message:
        "Retrieve reads for HLA genotyping by converting the alignments ({wildcards.nartype}) of group:{wildcards.group} from sample:{wildcards.sample} to reads in FASTQ format"
    shell:
        """
        mkdir -p tmp/
        (samtools sort -@ {threads} -n {input.reads} -T tmp/sort_{wildcards.sample}_{wildcards.group}_{wildcards.nartype}_ \
            | samtools fastq -@ {threads} >{output}) 2>{log}
        """
