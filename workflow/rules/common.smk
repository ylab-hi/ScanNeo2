import os

rawreads = {}
#samples = os.listdir(config["rnaseq"])
samples = [s for s in os.listdir(config['rnaseq']) if s.endswith('.'+config['rnaseq-format'])]
for i,v in enumerate(samples):
    rawreads["rep" + str(i+1)] = config["rnaseq"] + "/" + v

def get_fqs(wildcards):
    if config["rnaseq-format"] == "fastq":
        return rawreads[wildcards.sample]
    elif config["rnaseq-format"] == "bam":
        fqfile = "results/preproc/pre/"
        fqfile += wildcards.sample + ".fq.gz"
        print(fqfile)
        return fqfile

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
 


