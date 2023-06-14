import os
import re
from pathlib import Path

# check if files are a valid paired-end pair
def valid_paired_end(file1, file2):
    # check if both files are in FASTA format
    if Path(file1).suffix not in ['.fastq', '.fq'] and Path(file2).suffix not in ['.fastq', '.fq']:
        return False

    # check if first file contains _R1 or _fwd
    if not ("_R1" in file1 or "_fwd" in file1):
        return False
    # check if second file contains _R2 or _rev
    if not ("_R2" in file2 or "_rev" in file2):
        return False

    # check if the substrings until occurrence of either _R1/2 or _fwd/rev are equal
    file1_indicator = "_R1" if "_R1" in file1 else "_fwd"  
    file2_indicator = "_R2" if "_R2" in file2 else "_rev"
    file1_idx = file1.find(file1_indicator)
    file2_idx = file2.find(file2_indicator)

    if file1_idx != file2_idx:
        return False
    
    if file1[:file1_idx] != file2[:file2_idx]:
        return False

    # check if first file contains _R1 and the secon file contaif "_R1" in string1 and "_R2" not in string2:
        return False
    
    # check if first file contains _fwd and the secon file contains _rev
    if "_fwd" in string1 and "_rev" not in string2:
        return False

    return True


rnaseq_input = {}
rnaseq_filetype = None
rnaseq_files = config['rnaseq'].split(' ')
if len(rnaseq_files) == 1:
    if Path(rnaseq_files[0]).suffix in ['.fq', '.fastq', '.bam']:
        rnaseq_filetype = Path(rnaseq_files[0]).suffix # store file extension
        rnaseq_input[Path(rnaseq_files[0]).stem] = rnaseq_files[0]
    else:
        print("no rnaseq files found")
elif len(rnaseq_files) == 2:
    # check if both files are in valid paired-end format
    if valid_paired_end(rnaseq_files[0], rnaseq_files[1]):
        rnaseq_filetype = Path(rnaseq_files[0]).suffix # store file extension
        rnaseq_input.append((rnaseq_files[0], rnaseq_files[1]))
    else:
        print("no valid paired-end files found")

# determine input for star aligner
def star_output(wildcards):
    bamfiles = os.listdir("results/"+wildcards.sample+"/rnaseq/preproc/align/bams/")
    return bamfiles


# aggregate results from STAR alignment
def aggregate_star_align(wildcards):
    split_fastq_output = checkpoints.split_fastq.get(**wildcards).output[0]
    print(split_fastq_output)
    return expand("results/{sample}/rnaseq/preproc/align/bams/aln_{i}.bam",
        sample=wildcards.sample,
        i=glob_wildcards(os.path.join(split_fastq_output, "inputreads_{i}.fq.gz")).i)
    

def get_bams(wildcards):
    f = rnaseq_input[wildcards.sample]

    #if config["rnaseq-format"] == "fastq":
        #f = "results/preproc/"
        #f += wildcards.sample + "/"
        #f += wildcards.sample + ".bam"
    #elif config["rnaseq-format"] == "bam":
        #f = rawreads[wildcards.sample]

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
 



# find all samples in the rnaseq directory
#rnaseq_files = os.listdir(os.getcwd() + '/' + config['rnaseq'])
#print(config['rnaseq'])
#print(rnaseq_files)


#sample_counter = 1
#for i,v in enumerate(rnaseq_files):
    #if v.endswith('.bam'):
        #rnaseq_bam[v.split('.bam')[0]] = config['rnaseq'] + '/' + v
    #elif v.endswith(tuple(['.fastq','fq'])):
        #print("only fastq to be tested")
        ## check if paired-end reads
#        re = regex.compile(r'(.+)_R[12].fastq')
#        m = re.match(v)
#
#

