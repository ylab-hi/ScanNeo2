### align reads to genome using STAR (when reads are in FASTQ)
# Rule definitions are unconditional (issue #93): both this rule and the
# BAM-input chain (split_bamfile_RG -> ... -> merge_alignment_results) are
# always parsed, and they emit into filetype-tagged subdirectories (fq/ vs
# bam/) to avoid an AmbiguousRuleException. rnaseq_postproc_fixmate reads
# through get_rnaseq_star_bam to dispatch on the sample's filetype.
rule star_align_fastq:
    input:
        unpack(get_star_input),
        faidx="resources/refs/genome.fasta.fai",
        idx="resources/refs/star/",
    output:
        aln=temp("results/{sample}/rnaseq/align/fq/{group}_aligned_STAR.bam"),
        log="results/{sample}/rnaseq/align/fq/{group}_aligned_STAR.log",
        sj="results/{sample}/rnaseq/align/fq/{group}_aligned_STAR.tab",
    log:
        "logs/{sample}/align/star_align_fastq_{group}.log",
    threads: config["threads"]
    params:
        extra=lambda wildcards: (
            "--outSAMtype BAM Unsorted "
            "--genomeSAindexNbases 10 "
            "--outSAMattributes RG HI "
            f"--outSAMattrRGline ID:{wildcards.group} "
            "--outFilterMultimapNmax 50 "
            "--peOverlapNbasesMin 15 "
            "--alignSplicedMateMapLminOverLmate 0.5 "
            "--alignSJstitchMismatchNmax 5 -1 5 5 "
            "--chimOutType WithinBAM HardClip "
            f"--chimSegmentMin {config['align']['chimSegmentMin']} "
            f"--chimJunctionOverhangMin {config['align']['chimJunctionOverhangMin']} "
            f"--chimScoreDropMax {config['align']['chimScoreDropMax']} "
            f"--chimScoreMin {config['align']['chimScoreMin']} "
            "--chimScoreJunctionNonGTAG 0 "
            f"--chimScoreSeparation {config['align']['chimScoreSeparation']} "
            "--chimSegmentReadGapMax 3 "
            "--chimMultimapNmax 50 "
            "--outSAMstrandField intronMotif"
        ),
    message:
        "Aligning reads from {wildcards.group} to genome using STAR"
    wrapper:
        "v2.2.1/bio/star/align"


### align reads to genome using STAR (when reads are in BAM - no preprocessing performed)
checkpoint split_bamfile_RG:
    input:
        unpack(get_star_input),
    output:
        directory("results/{sample}/rnaseq/reads/{group}/bam/"),
    log:
        "logs/{sample}/align/split_bamfile_RG_{group}.log",
    conda:
        "../envs/samtools.yml"
    # samtools split is I/O-bound; per-thread gain plateaus past ~10.
    threads: min(10, config["threads"])
    shell:
        """
        mkdir -p {output}
        samtools split -@ {threads} \
            -u {output}/noRG.bam \
            -h {input} -f {output}/%!.%. {input} >{log} 2>&1
        """


rule bamfile_RG_to_fastq:
    input:
        "results/{sample}/rnaseq/reads/{group}/bam/{rg}.bam",
    output:
        "results/{sample}/rnaseq/reads/{group}/fastq/{rg}.fastq",
    log:
        "logs/{sample}/align/bamfile_RG_to_fastq_{group}_{rg}.log",
    conda:
        "../envs/samtools.yml"
    threads: config["threads"]
    message:
        "Converting group:{wildcards.group} BAM file to FASTQ for readgroup:{wildcards.rg}"
    shell:
        """
        (samtools collate -Oun128 -@ {threads} {input} \
            | samtools fastq -OT RG -@ {threads} - \
            | gzip -c - >{output}) 2>{log}
        """


rule star_align_bamfile:
    input:
        fq1="results/{sample}/rnaseq/reads/{group}/fastq/{rg}.fastq",
        faidx="resources/refs/genome.fasta.fai",
        idx="resources/refs/star/",
    output:
        aln="results/{sample}/rnaseq/align/{group}/{rg}.bam",
        log="results/{sample}/rnaseq/align/{group}/{rg}.log",
        sj="results/{sample}/rnaseq/align/{group}/{rg}.tab",
    log:
        "logs/{sample}/align/star_align_bamfile_{group}_{rg}.log",
    threads: config["threads"]
    params:
        extra=lambda wildcards: (
            "--outSAMtype BAM Unsorted --genomeSAindexNbases 10 "
            "--readFilesCommand zcat "
            f"--outSAMattributes RG HI --outSAMattrRGline ID:{wildcards.rg} "
            "--outFilterMultimapNmax 50 "
            "--peOverlapNbasesMin 15 "
            "--alignSplicedMateMapLminOverLmate 0.5 "
            "--alignSJstitchMismatchNmax 5 -1 5 5 "
            "--chimOutType WithinBAM HardClip "
            f"--chimSegmentMin {config['align']['chimSegmentMin']} "
            f"--chimJunctionOverhangMin {config['align']['chimJunctionOverhangMin']} "
            f"--chimScoreDropMax {config['align']['chimScoreDropMax']} "
            f"--chimScoreMin {config['align']['chimScoreMin']} "
            "--chimScoreJunctionNonGTAG 0 "
            f"--chimScoreSeparation {config['align']['chimScoreSeparation']} "
            "--chimSegmentReadGapMax 3 --chimMultimapNmax 50 "
            "--outSAMstrandField intronMotif"
        ),
    wrapper:
        "v2.2.1/bio/star/align"


rule merge_alignment_results:
    input:
        aggregate_aligned_rg,
    output:
        # `bam/` sub-path differentiates this from star_align_fastq so both
        # can always be defined; get_rnaseq_star_bam dispatches per sample
        # (issue #93).
        temp("results/{sample}/rnaseq/align/bam/{group}_aligned_STAR.bam"),
    log:
        "logs/{sample}/align/merge_alignment_results_{group}.log",
    threads: config["threads"]
    params:
        extra="",  # optional additional parameters as string
    wrapper:
        "v2.3.0/bio/samtools/merge"


# post-processingn and realignment
rule rnaseq_postproc_fixmate:
    input:
        get_rnaseq_star_bam,
    output:
        temp("results/{sample}/rnaseq/align/{group}_fixmate_STAR.bam"),
    log:
        "logs/{sample}/align/rnaseq_postproc_fixmate_{group}.log",
    conda:
        "../envs/samtools.yml"
    threads: 4
    params:
        mapq=f"--min-MQ {config['mapq']}",
    shell:
        """
        (
            samtools view -h -F 4 {params.mapq} {input} -o - \
                | samtools sort -n -@ {threads} -m4g -O SAM - -o - \
                | samtools fixmate -pcmu -O bam -@ {threads} - {output}
        ) >{log} 2>&1
        """


# sort and markdup needed to be separated (ensure no core dump for whatever reason)
rule rnaseq_postproc_markdup:
    input:
        bam="results/{sample}/rnaseq/align/{group}_fixmate_STAR.bam",
    output:
        "results/{sample}/rnaseq/align/{group}_final_STAR.bam",
    log:
        "logs/{sample}/align/rnaseq_postproc_markdup_{group}.log",
    conda:
        "../envs/samtools.yml"
    threads: 4
    resources:
        mem_mb_per_cpu=4000,
    shell:
        """
        mkdir -p tmp/
        samtools sort -@ {threads} -m4G -O BAM -T tmp/sort_{wildcards.sample}_{wildcards.group}_ {input.bam} \
            -o tmp/rnaseq_fixmate_sorted_{wildcards.sample}_{wildcards.group}.bam >{log} 2>&1
        samtools markdup -r -@ {threads} tmp/rnaseq_fixmate_sorted_{wildcards.sample}_{wildcards.group}.bam \
            {output} >>{log} 2>&1
        rm tmp/rnaseq_fixmate_sorted_{wildcards.sample}_{wildcards.group}.bam
        """


rule postproc_bam_index:
    input:
        "results/{sample}/rnaseq/align/{group}_final_STAR.bam",
    output:
        "results/{sample}/rnaseq/align/{group}_final_STAR.bam.bai",
    log:
        "logs/{sample}/align/postproc_bam_index_{group}.log",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools index {input} >{log} 2>&1
        """


## retrieve readgroups from bam file
rule get_readgroups:
    input:
        get_readgroups_input,
    output:
        "results/{sample}/{seqtype}/reads/{group}_readgroups.txt",
    log:
        "logs/{sample}/align/get_readgroups_{seqtype}_{group}.log",
    conda:
        "../envs/basic.yml"
    shell:
        """
        python workflow/scripts/get_readgroups.py '{input}' \
            {output} >{log} 2>&1
        """


# BWA-realign RNAseq alignments. {seqtype} is constrained to rnaseq so it
# doesn't compete with realign_dnaseq_bam / dnaseq_postproc for the
# canonical dnaseq output path (issue #93). The wildcard itself is kept
# because get_readgroups_input reads it.
rule realign:
    input:
        bam=get_readgroups_input,
        rg="results/{sample}/{seqtype}/reads/{group}_readgroups.txt",
        idx=multiext("resources/refs/bwa/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam="results/{sample}/{seqtype}/align/{group}_final_BWA.bam",
    log:
        "logs/{sample}/align/realign_{seqtype}_{group}.log",
    conda:
        "../envs/basic.yml"
    threads: config["threads"]
    wildcard_constraints:
        seqtype="rnaseq",
    shell:
        """
        (
            samtools collate -Oun128 {input.bam} \
                | samtools fastq -OT RG -@ {threads} - \
                | bwa mem -pt{threads} -CH <(cat {input.rg}) resources/refs/bwa/genome - - \
                | samtools sort -@6 -m1g -o {output}
        ) >{log} 2>&1
        """


# DNA-from-BAM input is BWA-realigned directly from the raw BAM (the OLD
# `realign` rule handled this via its {seqtype} wildcard; split out here so
# its output can carry a filetype tag). Output goes under bam/; the
# dnaseq_final_BWA_stage rule symlinks it to the canonical path that
# downstream consumers expect (issue #93).
rule realign_dnaseq_bam:
    input:
        bam=lambda wc: str(SAMPLES[wc.sample]["dnaseq"][wc.group]),
        rg="results/{sample}/dnaseq/reads/{group}_readgroups.txt",
        idx=multiext("resources/refs/bwa/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam="results/{sample}/dnaseq/align/bam/{group}_final_BWA.bam",
    log:
        "logs/{sample}/align/realign_dnaseq_bam_{group}.log",
    conda:
        "../envs/basic.yml"
    threads: config["threads"]
    shell:
        """
        (
            samtools collate -Oun128 {input.bam} \
                | samtools fastq -OT RG -@ {threads} - \
                | bwa mem -pt{threads} -CH <(cat {input.rg}) resources/refs/bwa/genome - - \
                | samtools sort -@6 -m1g -o {output}
        ) >{log} 2>&1
        """


# Materialize the canonical dnaseq _final_BWA.bam path from the filetype-tagged
# producer output (fq/ from dnaseq_postproc, bam/ from realign_dnaseq_bam) so
# downstream consumers (germline.smk, indel.smk, samtools_index_BWA_final)
# don't need a dispatcher of their own (issue #93).
rule dnaseq_final_BWA_stage:
    input:
        bam=get_dnaseq_final_bam_tagged,
    output:
        bam="results/{sample}/dnaseq/align/{group}_final_BWA.bam",
    log:
        "logs/{sample}/align/dnaseq_final_BWA_stage_{group}.log",
    conda:
        "../envs/basic.yml"
    shell:
        "ln -sfr {input.bam} {output.bam} 2>{log}"


### align DNAseq reads to genome using BWA (FASTQ input). Rule is unconditional
### per issue #93; get_dna_align_input decides per-sample which upstream to wire.
rule bwa_align_dnaseq:
    input:
        reads=get_dna_align_input,
        idx=multiext(
            "resources/refs/bwa/genome", ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    output:
        "results/{sample}/dnaseq/align/{group}_aligned_BWA.bam",
    log:
        "logs/{sample}/align/bwa_align_dnaseq_{group}.log",
    conda:
        "../envs/basic.yml"
    threads: config["threads"]
    params:
        extra="",
    shell:
        """
        (
            bwa mem -t{threads} resources/refs/bwa/genome \
                -R '@RG\\tID:{wildcards.group}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA' \
                {input.reads} \
                | samtools sort -@ 6 -n -m1g - -o {output}
        ) >{log} 2>&1
        """


rule dnaseq_postproc:
    input:
        aln="results/{sample}/dnaseq/align/{group}_aligned_BWA.bam",
    output:
        # `fq/` sub-path differentiates this from realign_dnaseq_bam so both
        # can always be defined; dnaseq_final_BWA_stage symlinks to the
        # canonical path per sample (issue #93).
        bam="results/{sample}/dnaseq/align/fq/{group}_final_BWA.bam",
    log:
        "logs/{sample}/align/dnaseq_postproc_{group}.log",
    conda:
        "../envs/samtools.yml"
    threads: 4
    resources:
        mem_mb=20000,
    params:
        extra="",
    shell:
        """
        mkdir -p tmp/
        (
            samtools fixmate -pcmu -O bam -@ {threads} {input.aln} - \
                | samtools sort -@ {threads} -m1g -O bam -T tmp/sort_{wildcards.sample}_{wildcards.group}_ - -o - \
                | samtools markdup -r -@ {threads} - {output.bam}
        ) >{log} 2>&1
        """


rule samtools_index_BWA_final:
    input:
        "results/{sample}/{seqtype}/align/{group}_final_BWA.bam",
    output:
        "results/{sample}/{seqtype}/align/{group}_final_BWA.bam.bai",
    log:
        "logs/{sample}/align/samtools_index_BWA_final_{seqtype}_{group}.log",
    threads: 4  # This value - 1 will be sent to -@
    params:
        extra="",  # optional params string
    wrapper:
        "v2.3.0/bio/samtools/index"
