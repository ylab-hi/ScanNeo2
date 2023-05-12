import os

rule bwa_index:
    input:
        fasta="refs/genome.fasta"
    output:
        idx = multiext("refs/bwa/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.26.0/bio/bwa/index"

rule realign:
    input:
        idx = multiext("refs/bwa/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        bam = "results/preproc/post/{sample}_dedup.bam"
    output:
        "results/indel/realign/{sample}.bam",
    log:
        "logs/bwa_realign_{sample}"
    conda:
        "../envs/realign.yml",
    threads: config['threads']
    shell:
        """
        samtools collate -Oun128 {input.bam} \
        | samtools fastq -OT RG,BC - \
        | bwa mem -pt8 -CH <(samtools view -H {input.bam}|grep ^@RG) misc/refgen/bwa/genome - \
        | samtools sort -@{threads} -m4g - \
        | samtools addreplacerg -r '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}' -o {output} -
        """

rule realign_index:
    input:
        "results/indel/realign/{sample}.bam",
    output:
        "results/indel/realign/{sample}.bam.bai",
    log:
        "logs/samtools_index_{sample}.log"
    params:
        extra="",  # optional params string
    threads: config['threads']  # This value - 1 will be sent to -@
    wrapper:
        "v1.26.0/bio/samtools/index"
        

rule transindel_build:
    input:
        bam = "results/indel/realign/{sample}.bam",
        idx = "results/indel/realign/{sample}.bam.bai"
    output:
        "results/indel/transindel/{sample}_redef.bam"
    log:
        "logs/transindel_buildbam_{sample}"
    conda:
        "../envs/transindel.yml"
    shell:
        """
        python3 workflow/scripts/transIndel/transIndel_build_RNA.py \
        -i {input.bam} \
        -o {output} \
        -r refs/genome.fasta \
        -g refs/genome.gtf > {log}
        """

rule transindel_call:
    input:
        "results/indel/transindel/{sample}_redef.bam"
    output:
        "results/indel/transindel/{sample}.indel.vcf"
    log:
        "logs/transindel_call_{sample}.log"
    conda:
        "../envs/transindel.yml"
    params:
        "-m {config[mapq]}"
    shell:
        """
        python workflow/scripts/transIndel/transIndel_call.py \
        -i {input} \
        -o results/indel/transindel/{wildcards.sample} \
        {params}
        """

rule slippage_removal:
    input:
        "results/indel/transindel/{sample}.indel.vcf"
    output:
        "results/indel/transindel/{sample}_slip.indel.vcf"
    conda:
        "../envs/transindel.yml"
    shell:
        """
        cp {input} {output}
        """

rule combine_indels:
    input:
        expand("results/indel/transindel/{sample}_slip.indel.vcf", sample=rawreads.keys())
    output:
        "results/indel/transindel/all.indel.vcf"
    shell:
        """
        cat {input} > {output}
        """

rule picard_create_dict:
    input:
        refgen =  config["refgen"],
    output:
        "misc/refgen/genome.dict"
    log:
        "logs/picard/create_dict.log",
    conda:
        "../envs/picard.yml"
    params:
        prefix=lambda wildcards, input: os.path.splitext(input[0])[0]
    shell:
        """
            if [ -f "{params}.dict" ]; then
                rm {params}.dict
            fi
            picard CreateSequenceDictionary -R {input} > {log}
            touch {output}
        """

# prepare file for haplotypecaller
rule gatk_realigner_target_creator:
    input:
        bam="results/indel/realign/{sample}.bam",
        bai="results/indel/realign/{sample}.bam.bai",
        ref=config["refgen"],
        fai="misc/refgen/genome.fa.fai",
        dict="misc/refgen/genome.dict",
#        known="dbsnp.vcf.gz",
#        known_idx="dbsnp.vcf.gz.tbi",
    output:
        intervals="results/indel/haplotypecaller/{sample}.intervals",
    log:
        "logs/gatk/realignertargetcreator/{sample}.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar -L misc/refgen/genome.dict",  # optional
    resources:
        mem_mb=1024,
    threads: 16
    wrapper:
        "v1.28.0/bio/gatk3/realignertargetcreator"


rule gatk_indel_realigner:
    input:
        bam="results/indel/realign/{sample}.bam",
        bai="results/indel/realign/{sample}.bam.bai",
        ref=config["refgen"],
        fai="misc/refgen/genome.fa.fai",
        dict="misc/refgen/genome.dict",
        target_intervals="results/indel/haplotypecaller/{sample}.intervals",
    output:
        bam="results/indel/haplotypecaller/{sample}.realigned.bam",
        bai="results/indel/haplotypecaller/{sample}.realigned.bai",
    log:
        "logs/gatk3/indelrealigner/{sample}.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    threads: 16
    resources:
        mem_mb=1024,
    wrapper:
        "v1.28.0/bio/gatk3/indelrealigner"


rule gatk_haplotype_caller_gvcf:
    input:
#        bam="results/indel/haplotypecaller/{sample}.realigned.bam",
        bam="results/preproc/post/{sample}_dedup.bam",
        ref=config["refgen"]
    output:
        gvcf="results/indel/haplotypecaller/{sample}.g.vcf",
        bam="results/indel/haplotypecaller/{sample}.assemb_haplo.bam",
    log:
        "logs/gatk/haplotypecaller/{sample}.log",
    params:
        extra="--sequence-dictionary misc/refgen/genome.dict",  # optional
        java_opts="",  # optional
    threads: config['threads']
    resources:
        mem_mb=1024,
    wrapper:
        "v1.27.0/bio/gatk/haplotypecaller"

rule combine_gvcfs:
    input:
        gvcfs=expand("results/indel/haplotypecaller/{sample}.g.vcf", sample=config['rnaseq']),
        ref=config["refgen"],
    output:
        gvcf="results/indel/haplotypecaller/all.g.vcf",
    log:
        "logs/gatk/combinegvcfs.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.28.0/bio/gatk/combinegvcfs"

rule genotype_gvcfs:
    input:
        gvcf="results/indel/haplotypecaller/all.g.vcf",
    # N.B. gvcf or genomicsdb must be specified
    # in the latter case, this is a GenomicsDB data store
        ref=config["refgen"]
    output:
        vcf="results/indel/haplotypecaller/all.vcf",
    log:
        "logs/gatk/genotypegvcfs.log"
    params:
        extra="",  # optional
        java_opts="", # optional
    resources:
        mem_mb=1024
    wrapper:
        "v1.28.0/bio/gatk/genotypegvcfs"


rule picard_split_vcfs:
    input:
        "results/indel/haplotypecaller/all.vcf"
    output:
        snp="results/indel/snp.vcf",
        indel="results/haplotypecaller/indel.vcf"
    log:
        "logs/spltvcfs.log"
    conda:
        "../envs/picard.yml"
    shell:
        """
            picard SplitVcfs I={input} \
            SNP_OUTPUT={output.snp} \
            INDEL_OUTPUT={output.indel}
            STRICT=false
        """


rule merge_indels:
    input:
        #indel="results/haplotypecaller/indel.vcf",
        transindel="results/indel/transindel/all.indel.vcf"
    output:
        "results/indel/indel.vcf"
    shell:
        """
        cp {input} {output}
        """

