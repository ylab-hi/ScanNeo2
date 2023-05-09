def comma_join(file_list):
    return ",".join(file_list)


rule spladder:
    input:
        bam = "results/preproc/post/{sample}_dedup.bam",
        bamidx = "results/preproc/post/{sample}_dedup.bam.bai",
    output:
        "results/as/{sample}/merge_graphs.txt"
    conda: 
        "../envs/spladder.yml"
    log:
        "logs/as/spladder_{sample}.log"
    params:
        bams = comma_join(expand("results/preproc/post/{sample}_dedup.bam", sample=config["rnaseq"])),
    shell:
        """
        spladder build -b {params.bams} \
        -a {config[annotation]} -o {output} \
        --filter-overlap-exons --no-primary-only \
        --confidence 2 --quantify-graph \
        --qmode all 2> {log}
        """

rule splicing_to_vcf:
    input:
        "results/as/{sample}/merge_graphs.txt"
    output:
        "results/as/{sample}.vcf"
    conda:
        "../envs/spladder.yml"
    shell:
        """
        """

rule combine_splicing:
    input:
        expand("results/as/{sample}.vcf", sample=config["rnaseq"])
    output:
        "results/as/all.vcf"
    shell:
        """
        """


#TYPE = ["alt_3prime", "alt_5prime", "exon_skip", "intron_retention", "mutex_exons"] 
#TYPETT = ["A3","A5","ES","IR","MUT"]
#type_conv = {
#    "A3": "alt_3prime",
#    "A5": "alt_5prime",
#    "ES": "exon_skip",
#    "IR": "intron_retention",
#    "MUT": "mutex_exons"
#}

#type_to_param = {
#    "alt_3prime": "A3", 
#    "alt_5prime": "A5", 
#    "exon_skip": "ES", 
#    "intron_retention": "IR", 
#    "mutex_exons": "MUT"}

#rule bisbee_counts:
#    input:
#        lambda wildcards: "results/splicing/spladder/{}/merge_graphs_{}_C2.counts.hdf5".format(wildcards.sample, type_conv[wildcards.type])
#    output:
#        "results/splicing/bisbee/{sample}/counts/{sample}.{type}.bisbeeCounts.csv"
#    conda: 
#        "../envs/bisbee.yml"
#    log:
#        "logs/bisbee_prep_{sample}_{type}.log"
#    shell:
#        "python3 workflow/scripts/bisbee/utils/prep.py {input} {wildcards.type} results/splicing/bisbee/{wildcards.sample}/counts/{wildcards.sample}.{wildcards.type} 3 2> {log}"


#rule bisbee_outliers:
#    input:
#        "results/splicing/bisbee/{sample}/counts/{sample}.{type}.bisbeeCounts.csv"
#    output:
#        fit = "results/splicing/bisbee/{sample}/outlier/{sample}.{type}.bisbeeFit.csv",
#        score = "results/splicing/bisbee/{sample}/outlier/{sample}.{type}.bisbeeOutlier.csv"
#    message:
#        "commence outlier analysis"
#    conda:
#        "../envs/bisbee_outliers.yml"
#    log:
#        "logs/bisbee_outliersFit_{sample}_{type}.log",
#        "logs/bisbee_outliersScore_{sample}_{type}.log"
#    shell:
#        """
#        Rscript workflow/scripts/bisbee/stats/outlierFit.R \
#        {input} 80 \
#        results/splicing/bisbee/{wildcards.sample}/outlier/{wildcards.sample}.{wildcards.type} \
#        > {log[0]}
#        Rscript workflow/scripts/bisbee/stats/outlierScore.R \
#        {output.fit} {input} \
#        results/splicing/bisbee/{wildcards.sample}/outlier/{wildcards.sample}.{wildcards.type} \
#        > {log[1]}
#        """

#rule bisbee_outliers_filter:
#    input:
#        "results/splicing/bisbee/{sample}/outlier/{sample}.{type}.bisbeeOutlier.csv"
#    output:
#        main = "results/splicing/bisbee/{sample}/outlier/{sample}.{type}.bisbeeOutlier.thresh10.min1test.csv",
#        scores =  "results/splicing/bisbee/{sample}/outlier/{sample}.{type}.bisbeeOutlier.thresh10.min1test.scores.csv",
#        summary = "results/splicing/bisbee/{sample}/outlier/{sample}.{type}.bisbeeOutlier.thresh10.min1test.summary.csv"
#    conda: 
#        "../envs/bisbee.yml"
#    log:
#        "logs/bisbee_outlierfilter_{sample}_{type}.log"
#    shell:
#        """
#        python3 workflow/scripts/bisbee/utils/filtOut.py \
#        results/splicing/bisbee/{wildcards.sample}/outlier/ \
#        {wildcards.sample}.{wildcards.type} 10 1
#        """
