import os

rawreads = {}
#samples = os.listdir(config["rnaseq"])
samples = [s for s in os.listdir() if s.endswith('.'+config['rnaseq-format'])]
for i,v in enumerate(samples):
    rawreads["rep" + str(i+1)] = config["rnaseq"] + "/" + v

def get_fqs(wildcards):
    return rawreads[wildcards.sample]

def get_bams(wildcards):
    if config["rnaseq-format"] == "fastq":
        f = "results/preproc/"
        f += wildcards.sample + "/"
        f += wildcards.sample + ".bam"
    elif config["rnaseq-format"] == "bam":
        f = rawreads[wildcards.sample]

    return f

def change_ext(filename, ext):
    return os.path.splitext(filename)[0] + "." + ext

def file_exists(path):
    return os.path.exists(path)

rule combine_sources:
    input:
        "results/indel/transindel/all.indel.vcf",
    output:
        "results/variants/variants.vcf"
    shell:
        "cat {input} > {output}"

rule variants_gzip:
    input:
        "results/variants/variants.vcf"
    output:
        "results/variants/variants.vcf.gz",
    log:
        "logs/variants/index.log"
    conda:
        "../envs/samtools.yml"
    shell:
        "bgzip -c {input} > {output}"


rule download_vep_plugins:
    output:
        directory("resources/vep/plugins")
    params:
        release=100
    wrapper:
        "v1.29.0/bio/vep/plugins"

rule variant_effect_predictor:
    input: 
        calls="results/variants/variants.vcf.gz",
        plugins="resources/vep/plugins/",
        fasta="refs/genome.fasta",
        fai="refs/genome.fasta.fai",
        gff="refs/genome.gtf.gz",
        csi="refs/genome.gtf.gz.csi"
    output:
        calls="results/variants/variants.annotated.vcf",  # .vcf, .vcf.gz or .bcf
        stats="results/variants/variants.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["NMD"],
        extra="--everything",  # optional: extra arguments
    log:
        "logs/vep/annotate.log",
    threads: 8
    wrapper:
        "v1.29.0/bio/vep/annotate"
    
        




#def get_refdict(wildcards):
    #refgen = config["refgen"]
    #base = os.path.splitext(refgen)[0]
    #return os.rename(refgen, base + '.dict')


#def num_samples(wildcards):
#    return len(config["samples"][wildcards.sample])


# rule to merge VCFs


#rule merge_vcf:
#    input:
#        "results/exitron/final.vcf"
#        "results/splicing/final.vcf"
#        "results/genefusion/final.vcf"
#    output:
#        "results/variants.vcf"
#
#    shell:
#        "workflow/scripts/merge_vcf.py {config[refgen]} {input} {output} 2> logs/error.err"
#        
 


