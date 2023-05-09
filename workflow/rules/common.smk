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


rule variant_effect_predictor:
    input: 
        "results/genefusion/all.vcf", 
        "results/as/all.vcf", 
        "results/exitron/all.vcf", 
        "results/indel/indel.vcf",
        "results/indel/snp.vcf"
    output: "results/vep/output.tab"
    shell:
        """
        """
















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
 


