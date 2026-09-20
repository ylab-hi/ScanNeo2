import os
from snakemake.remote import HTTP


# transindel BAM-rebuild streams the whole BAM single-threaded (~2 h on deep RNA),
# so it is scattered by chromosome: split the aligned BAM, rebuild each chromosome
# in parallel, then merge. Each read's redefined CIGAR depends only on that read,
# the reference, and the genome-wide GTF splice coverage -- no cross-read state --
# so per-chromosome rebuild + merge is equivalent to a whole-BAM rebuild.
checkpoint split_bam_ti_build:
    input:
        bam="results/{sample}/{seqtype}/align/{group}_final_BWA.bam",
        idx="results/{sample}/{seqtype}/align/{group}_final_BWA.bam.bai",
    output:
        temp(
            directory(
                "results/{sample}/{seqtype}/indel/transindel/{group}_build_split"
            )
        ),
    log:
        "logs/{sample}/indel/ti_split_{seqtype}_{group}.log",
    conda:
        "../envs/basic.yml"
    message:
        "Splitting BAM by chromosome for transindel build on sample:{wildcards.sample} with group:{wildcards.group}"
    shell:
        """
        python workflow/scripts/split_bam_by_chr.py {input.bam} {output} >{log} 2>&1
        """


rule detect_long_indel_ti_build_RNA:
    input:
        bam="results/{sample}/rnaseq/indel/transindel/{group}_build_split/{chr}.bam",
    output:
        bam=temp(
            "results/{sample}/rnaseq/indel/transindel/{group}_build_perchr/{chr}.bam"
        ),
    log:
        "logs/{sample}/indel/ti_build_RNA_{group}_{chr}.log",
    conda:
        "../envs/transindel.yml"
    message:
        "transindel build (RNA) on sample:{wildcards.sample} group:{wildcards.group} chr:{wildcards.chr}"
    shell:
        """
        python3 workflow/scripts/transindel/transIndel_build_RNA.py \
            -i {input.bam} \
            -o {output.bam} \
            -r resources/refs/genome.fasta \
            -g resources/refs/genome.gtf >{log} 2>&1
        """


rule detect_long_indel_ti_build_DNA:
    input:
        bam="results/{sample}/dnaseq/indel/transindel/{group}_build_split/{chr}.bam",
    output:
        bam=temp(
            "results/{sample}/dnaseq/indel/transindel/{group}_build_perchr/{chr}.bam"
        ),
    log:
        "logs/{sample}/indel/ti_build_DNA_{group}_{chr}.log",
    conda:
        "../envs/transindel.yml"
    message:
        "transindel build (DNA) on sample:{wildcards.sample} group:{wildcards.group} chr:{wildcards.chr}"
    shell:
        """
        python workflow/scripts/transindel/transIndel_build_DNA.py \
            -i {input.bam} -o {output.bam} >{log} 2>&1
        """


rule detect_long_indel_ti_call:
    input:
        bam="results/{sample}/{seqtype}/indel/transindel/{group}_build_perchr/{chr}.bam",
    output:
        temp(
            "results/{sample}/{seqtype}/indel/transindel/{group}_call_perchr/{chr}.indel.vcf"
        ),
    log:
        "logs/{sample}/indel/ti_call_{seqtype}_{group}_{chr}.log",
    conda:
        "../envs/transindel.yml"
    params:
        mapq=config["mapq"],
    # transIndel_call scans the BAM sequentially via pileup(region=None), so the
    # per-chromosome split needs no index. It is called per-chr rather than once
    # over the whole build BAM so the pileup parallelizes across the cluster.
    message:
        "Calling long indels with transindel on sample:{wildcards.sample} group:{wildcards.group} chr:{wildcards.chr}"
    shell:
        """
        python workflow/scripts/transindel/transIndel_call.py \
            -i {input.bam} \
            -l 10 \
            -o results/{wildcards.sample}/{wildcards.seqtype}/indel/transindel/{wildcards.group}_call_perchr/{wildcards.chr} \
            -m {params} >{log} 2>&1
        """


rule merge_ti_call:
    input:
        aggregate_ti_call,
    output:
        "results/{sample}/{seqtype}/indel/transindel/{group}_call.indel.vcf",
    log:
        "logs/{sample}/indel/ti_merge_call_{seqtype}_{group}.log",
    conda:
        "../envs/basic.yml"
    message:
        "Merging per-chromosome transindel calls on sample:{wildcards.sample} with group:{wildcards.group}"
    # Concatenate the per-chr transindel VCFs into one call set under a single
    # canonical sample column (transIndel names each per-file sample after its
    # -o prefix, so they otherwise disagree). This matches a whole-BAM
    # transIndel_call: contig lines and the final coordinate sort are added
    # downstream by long_indel_augment / longindel_sort_and_compress.
    shell:
        """
        (
            first=$(echo {input} | tr ' ' '\\n' | head -n1)
            grep '^##' "$first"
            printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tresults/{wildcards.sample}/{wildcards.seqtype}/indel/transindel/{wildcards.group}_call\\n'
            grep -hv '^#' {input} || [ $? -eq 1 ]
        ) >{output} 2>{log}
        """


# resove alleles and remove PCR slippage
rule long_indel_augment:
    input:
        "results/{sample}/{seqtype}/indel/transindel/{group}_call.indel.vcf",
    output:
        "results/{sample}/{seqtype}/indel/transindel/{group}_long.indels_augmented.vcf",
    log:
        "logs/{sample}/indel/long_indel_augment_{seqtype}_{group}.log",
    conda:
        "../envs/manipulate_vcf.yml"
    message:
        "Augment long indels with group and source information and resolving alleles and removing PCR slippage using transindel on sample:{wildcards.sample} group:{wildcards.group}"
    shell:
        """
        python3 workflow/scripts/add_infos_to_vcf.py \
            {input} \
            long_indel \
            {wildcards.group} \
            {output}_infos >{log} 2>&1

        python3 workflow/scripts/add_contigs_to_vcf.py \
            {output}_infos \
            {output}_contigs \
            resources/refs/genome.fasta >>{log} 2>&1

        python3 workflow/scripts/slippage_removal.py \
            resources/refs/genome.fasta \
            {output}_contigs \
            {output} >>{log} 2>&1

        rm -r {output}_contigs {output}_infos
        """


rule longindel_sort_and_compress:
    input:
        "results/{sample}/{seqtype}/indel/transindel/{group}_long.indels_augmented.vcf",
    output:
        "results/{sample}/{seqtype}/indel/transindel/{group}_long.indels.vcf.gz",
    log:
        "logs/{sample}/indel/longindel_sort_{seqtype}_{group}.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Sorting and compressing long indels on sample:{wildcards.sample}"
    shell:
        """
        (bcftools sort {input} -o - | bcftools view -O z -o {output}) >{log} 2>&1
        """


rule combine_longindels:
    input:
        get_longindels,
    output:
        "results/{sample}/variants/long.indels.vcf.gz",
    log:
        "logs/{sample}/indel/combine_longindels.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Combining long indels on sample:{wildcards.sample}"
    shell:
        """
        (bcftools concat --naive-force -O z {input} -o - | bcftools sort -O z -o {output}) >{log} 2>&1
        """


####### MUTECT2 ######


checkpoint split_bam_detect_short_indels_m2:
    input:
        bam="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.baserecal.bam",
        idx="results/{sample}/{seqtype}/indel/htcaller/{group}_variants.1rd.baserecal.bam.bai",
    output:
        directory("results/{sample}/{seqtype}/indel/mutect2/{group}_baserecal_split"),
    log:
        "logs/{sample}/indel/split_bam_m2_{seqtype}_{group}.log",
    conda:
        "../envs/basic.yml"
    message:
        "Splitting bam file for somatic SNV/Indel detection with Mutect2 on recalibrated data on sample:{wildcards.sample} with group:{wildcards.group}"
    shell:
        """
        python workflow/scripts/split_bam_by_chr.py \
            {input.bam} {output} >{log} 2>&1
        """


rule detect_short_indels_m2:
    input:
        # the .bai sits next to this .bam inside the split checkpoint directory
        # (written by split_bam_by_chr.py); GATK finds it automatically
        map="results/{sample}/{seqtype}/indel/mutect2/{group}_baserecal_split/{chr}.bam",
        # matched-normal BAM (+index) when present -> paired Mutect2; [] otherwise.
        # The wrapper adds it to the command via params.extra (`-I ... -normal`).
        normal=get_mutect_normal_input,
        fasta="resources/refs/genome.fasta",
    output:
        vcf=temp(
            "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/raw/{chr}.vcf"
        ),
        # Mutect2 writes this sidecar next to the vcf; FilterMutectCalls requires
        # it. Declare it explicitly so snakemake tracks/stages it alongside the
        # vcf (the raw/ subdir split it from filter's working directory).
        stats=temp(
            "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/raw/{chr}.vcf.stats"
        ),
        # per-chr read-orientation counts; aggregated by learn_read_orientation_m2
        # into the group's artifact-priors model for FilterMutectCalls
        f1r2=temp(
            "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/raw/{chr}.f1r2.tar.gz"
        ),
    log:
        "logs/{sample}/indel/detect_short_indels_m2_{seqtype}_{group}_{chr}.log",
    threads: 4
    resources:
        mem_mb=1024,
    params:
        extra=get_mutect_paired_extra,
    message:
        "Detection of somatic SNVs/Indels with Mutect2 on sample:{wildcards.sample} with group:{wildcards.group} on chromosome {wildcards.chr}"
    wrapper:
        "v1.31.1/bio/gatk/mutect"


rule learn_read_orientation_m2:
    input:
        f1r2=aggregate_f1r2_mutect2,
    output:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_read-orientation-model.tar.gz",
    log:
        "logs/{sample}/indel/learn_read_orientation_m2_{seqtype}_{group}.log",
    conda:
        "../envs/gatk.yml"
    resources:
        mem_mb=4096,
    message:
        "Learning read-orientation model (FFPE/OxoG artifact priors) on sample:{wildcards.sample} with group:{wildcards.group}"
    shell:
        # LearnReadOrientationModel takes one -I per per-chromosome f1r2 archive
        """
        (
            tmp=$(mktemp -d)
            trap 'st=$?; rm -rf "$tmp" || true; exit $st' EXIT
            gatk LearnReadOrientationModel $(printf -- '-I %s ' {input.f1r2}) \
                -O {output} --tmp-dir "$tmp"
        ) >{log} 2>&1
        """


rule filter_short_indels_m2:
    input:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants/raw/{chr}.vcf",
        stats="results/{sample}/{seqtype}/indel/mutect2/{group}_variants/raw/{chr}.vcf.stats",
        bam="results/{sample}/{seqtype}/indel/mutect2/{group}_baserecal_split/{chr}.bam",
        ref="resources/refs/genome.fasta",
        # the wrapper maps input.f1r2 -> --orientation-bias-artifact-priors
        f1r2="results/{sample}/{seqtype}/indel/mutect2/{group}_read-orientation-model.tar.gz",
    output:
        vcf=temp(
            "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf"
        ),
    log:
        "logs/{sample}/indel/filter_short_indels_m2_{seqtype}_{group}_{chr}.log",
    resources:
        mem_mb=1024,
    params:
        extra=(
            "--max-alt-allele-count 3 "
            f"--min-median-base-quality {config['basequal']} "
            f"--min-median-mapping-quality {config['mapq']} "
            f"--threshold-strategy {config['indel']['strategy']} "
            f"--f-score-beta {config['indel']['fscorebeta']} "
            f"--false-discovery-rate {config['indel']['fdr']} "
            f"--pcr-slippage-rate {config['indel']['sliprate']} "
            f"--min-slippage-length {config['indel']['sliplen']}"
        ),
        java_opts="",  # optional
    message:
        "Filtering somatic SNVs/Indels with FilterMutectCalls on sample:{wildcards.sample} with group:{wildcards.group} on chromosome {wildcards.chr}"
    wrapper:
        "v1.31.1/bio/gatk/filtermutectcalls"


rule sort_short_indels_m2:
    input:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf",
    output:
        temp(
            "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf.gz"
        ),
    log:
        "logs/{sample}/indel/sort_short_indels_m2_{seqtype}_{group}_{chr}.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Sorting vcf file from somatic variant calling (mutect2) on recalibrated data on sample:{wildcards.sample} with group:{wildcards.group} on chromosome {wildcards.chr}"
    shell:
        """
        (bcftools sort {input} -o - | bcftools view -O z -o {output}) >{log} 2>&1
        """


rule index_short_indels_m2:
    input:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf.gz",
    output:
        temp(
            "results/{sample}/{seqtype}/indel/mutect2/{group}_variants/{chr}_flt.vcf.gz.tbi"
        ),
    log:
        "logs/{sample}/indel/index_short_indels_m2_{seqtype}_{group}_{chr}.log",
    message:
        "Indexing vcf file from somatic variant valling (mutect2) first round on recalibrated data on sample:{wildcards.sample} with group:{wildcards.group} on chromosome {wildcards.chr}"
    wrapper:
        "v4.0.0/bio/bcftools/index"


rule merge_short_indels_m2:
    input:
        calls=aggregate_vcf_mutect2,
        idx=aggregate_idx_mutect2,
    output:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz",
    log:
        "logs/{sample}/indel/merge_short_indels_m2_{seqtype}_{group}.log",
    params:
        extra="-a",
    message:
        "Merging vcf files from first round of variant calling (htcaller) on original, unrecalibrated data on sample:{wildcards.sample} with group:{wildcards.group}"
    wrapper:
        "v4.0.0/bio/bcftools/concat"


######### POST-PROCESSING ########


rule index_merged_short_indels_m2:
    input:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz",
    output:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz.tbi",
    log:
        "logs/{sample}/indel/index_merged_short_indels_m2_{seqtype}_{group}.log",
    message:
        "Indexing vcf file from somatic variant valling (mutect2) on merged recalibrated data on sample:{wildcards.sample} with group:{wildcards.group}"
    wrapper:
        "v4.0.0/bio/bcftools/index"


rule select_short_indels_m2:
    input:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz",
        idx="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz.tbi",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf",
    log:
        "logs/{sample}/indel/select_short_indels_m2_{seqtype}_{group}.log",
    resources:
        mem_mb=1024,
    params:
        # --exclude-filtered: keep only PASS calls. FilterMutectCalls only
        # annotates the FILTER column, so without this the non-PASS calls
        # (weak_evidence, orientation, clustered_events, ...) flow downstream and
        # are treated as somatic neoepitope sources.
        extra="--select-type-to-include INDEL --exclude-filtered",
        java_opts="",  # optional
    message:
        "Selecting short somatic indels with SelectVariants on sample:{wildcards.sample}"
    wrapper:
        "v1.31.1/bio/gatk/selectvariants"


rule augment_short_indels_m2:
    input:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf",
    output:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels_augmented.vcf",
    log:
        "logs/{sample}/indel/augment_short_indels_m2_{seqtype}_{group}.log",
    conda:
        "../envs/manipulate_vcf.yml"
    message:
        "Combining somatic short indels detected by Mutect2 on sample:{wildcards.sample}"
    shell:
        """
        python workflow/scripts/add_infos_to_vcf.py \
            {input} \
            short_indel \
            {wildcards.group} \
            {output} >{log} 2>&1
        """


rule sort_aug_short_indels_m2:
    input:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels_augmented.vcf",
    output:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf.gz",
    log:
        "logs/{sample}/indel/sort_aug_short_indels_m2_{seqtype}_{group}.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Sorting and compressing short indels on sample:{wildcards.sample}"
    # Keep only the tumor sample. Paired dnaseq Mutect2 emits a matched-normal
    # column alongside the tumor, but combine_aug_short_indels_m2 stacks this
    # dnaseq VCF with the tumor-only rnaseq one, which requires both to carry the
    # same single tumor sample.
    shell:
        """
        (bcftools sort {input} -o - | bcftools view -s {wildcards.sample}_{wildcards.group} -O z -o {output}) >{log} 2>&1
        """


rule combine_aug_short_indels_m2:
    input:
        get_shortindels,
    output:
        "results/{sample}/variants/somatic.short.indels.vcf.gz",
    log:
        "logs/{sample}/indel/combine_short_indels.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Combining short indels on sample:{wildcards.sample}"
    shell:
        """
        (bcftools concat --naive-force -O z {input} -o - | bcftools sort -O z -o {output}) >{log} 2>&1
        """


rule select_SNVs_m2:
    input:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz",
        idx="results/{sample}/{seqtype}/indel/mutect2/{group}_variants.vcf.gz.tbi",
        ref="resources/refs/genome.fasta",
    output:
        vcf="results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf",
    log:
        "logs/{sample}/indel/select_SNVs_m2_{seqtype}_{group}.log",
    resources:
        mem_mb=1024,
    params:
        # --exclude-filtered: keep only PASS calls (see select_short_indels_m2)
        extra="--select-type-to-include SNP --exclude-filtered",
        java_opts="",  # optional
    message:
        "Selecting somatic SNVs with SelectVariants on sample:{wildcards.sample}"
    wrapper:
        "v1.31.1/bio/gatk/selectvariants"


rule augment_somatic_SNVs_m2:
    input:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf",
    output:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs_augmented.vcf",
    log:
        "logs/{sample}/indel/augment_somatic_SNVs_m2_{seqtype}_{group}.log",
    conda:
        "../envs/manipulate_vcf.yml"
    message:
        "Combining somatic SNVs detected by Mutect2 on sample:{wildcards.sample}"
    shell:
        """
        python workflow/scripts/add_infos_to_vcf.py \
            {input} \
            snv \
            {wildcards.group} \
            {output} >{log} 2>&1
        """


rule sort_somatic_SNVs_m2:
    input:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs_augmented.vcf",
    output:
        "results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf.gz",
    log:
        "logs/{sample}/indel/sort_somatic_SNVs_m2_{seqtype}_{group}.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Sorting and compressing somatic SNVs on sample:{wildcards.sample}"
    # Keep only the tumor sample (see sort_aug_short_indels_m2): combine_somatic_SNVs_m2
    # likewise stacks the paired dnaseq VCF with the tumor-only rnaseq one.
    shell:
        """
        (bcftools sort {input} -o - | bcftools view -s {wildcards.sample}_{wildcards.group} -O z -o {output}) >{log} 2>&1
        """


rule combine_somatic_SNVs_m2:
    input:
        get_snvs,
    output:
        "results/{sample}/variants/somatic.snvs.vcf.gz",
    log:
        "logs/{sample}/indel/combine_somatic_SNVs.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Combining somatic SNVs on sample:{wildcards.sample}"
    shell:
        """
        (bcftools concat --naive-force -O z {input} -o - | bcftools sort -O z -o {output}) >{log} 2>&1
        """


######### RNA GERMLINE SUBTRACTION ########
# RNA is called tumor-only, so its germline is not removed during calling. When
# the sample has a matched normal, subtract that normal's germline calls (the
# final-round HaplotypeCaller/VQSR set) from the RNA somatic calls. DNA is
# already germline-subtracted by paired Mutect2 and does not route through here.


rule build_germline_reference:
    input:
        snvs="results/{sample}/{seqtype}/indel/htcaller/{group}_snvs.final.flt.vcf",
        indels="results/{sample}/{seqtype}/indel/htcaller/{group}_indel.final.flt.vcf",
    output:
        vcf="results/{sample}/{seqtype}/indel/htcaller/{group}_germline.final.vcf.gz",
        tbi="results/{sample}/{seqtype}/indel/htcaller/{group}_germline.final.vcf.gz.tbi",
    log:
        "logs/{sample}/indel/build_germline_reference_{seqtype}_{group}.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Building germline reference (final-round SNVs + indels) for sample:{wildcards.sample} group:{wildcards.group}"
    shell:
        """
        (
            tmp=$(mktemp -d)
            bcftools view -O z -o $tmp/snvs.vcf.gz {input.snvs} && bcftools index -t $tmp/snvs.vcf.gz
            bcftools view -O z -o $tmp/indels.vcf.gz {input.indels} && bcftools index -t $tmp/indels.vcf.gz
            bcftools concat -a $tmp/snvs.vcf.gz $tmp/indels.vcf.gz | bcftools sort -O z -o {output.vcf}
            bcftools index -t {output.vcf}
            rm -rf $tmp
        ) >{log} 2>&1
        """


rule subtract_germline_snvs:
    input:
        vcf="results/{sample}/rnaseq/indel/mutect2/{group}_somatic.snvs.vcf.gz",
        germline=matched_normal_germline_vcf,
        germline_idx=matched_normal_germline_tbi,
    output:
        "results/{sample}/rnaseq/indel/mutect2/{group}_somatic.snvs.germsub.vcf.gz",
    log:
        "logs/{sample}/indel/subtract_germline_snvs_{group}.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Subtracting matched-normal germline from RNA somatic SNVs on sample:{wildcards.sample} group:{wildcards.group}"
    shell:
        """
        (
            tmp=$(mktemp -d)
            cp {input.vcf} $tmp/in.vcf.gz && bcftools index -t $tmp/in.vcf.gz
            bcftools isec -C -w1 -O z -o {output} $tmp/in.vcf.gz {input.germline}
            rm -rf $tmp
        ) >{log} 2>&1
        """


rule annotate_rna_editing:
    input:
        vcf="results/{sample}/rnaseq/indel/mutect2/{group}_somatic.snvs.germsub.vcf.gz",
        redi="resources/rediportal/rediportal_hg38.txt.gz",
        redi_idx="resources/rediportal/rediportal_hg38.txt.gz.tbi",
    output:
        "results/{sample}/rnaseq/indel/mutect2/{group}_somatic.snvs.germsub.reanno.vcf.gz",
    log:
        "logs/{sample}/indel/annotate_rna_editing_{group}.log",
    conda:
        "../envs/basic.yml"
    message:
        "Annotating RNA somatic SNVs against REDIportal (A-to-I editing) on sample:{wildcards.sample} group:{wildcards.group}"
    shell:
        """
        python workflow/scripts/annotate_rna_editing.py \
            {input.vcf} {input.redi} {output} >{log} 2>&1
        """


rule subtract_germline_short_indels:
    input:
        vcf="results/{sample}/rnaseq/indel/mutect2/{group}_somatic.short.indels.vcf.gz",
        germline=matched_normal_germline_vcf,
        germline_idx=matched_normal_germline_tbi,
    output:
        "results/{sample}/rnaseq/indel/mutect2/{group}_somatic.short.indels.germsub.vcf.gz",
    log:
        "logs/{sample}/indel/subtract_germline_short_indels_{group}.log",
    conda:
        "../envs/bcftools.yml"
    message:
        "Subtracting matched-normal germline from RNA somatic short indels on sample:{wildcards.sample} group:{wildcards.group}"
    shell:
        """
        (
            tmp=$(mktemp -d)
            cp {input.vcf} $tmp/in.vcf.gz && bcftools index -t $tmp/in.vcf.gz
            bcftools isec -C -w1 -O z -o {output} $tmp/in.vcf.gz {input.germline}
            rm -rf $tmp
        ) >{log} 2>&1
        """
