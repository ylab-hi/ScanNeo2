import os
import re
import glob
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


# aggregate results from STAR alignment
def aggregate_alignments(wildcards):
    # make sure that all samples are processed in checkpoint - split fastq file
    checkpoint_output = checkpoints.split_bamfile_RG.get(**wildcards).output[0]
    return expand("results/{sample}/rnaseq/align/bamfiles/{readgroup}.bam",
        sample=wildcards.sample,
        readgroup=glob_wildcards(os.path.join(checkpoint_output, "{readgroup}.bam")).readgroup)


def aggregate_bwa_alignments(wildcards):
    split_bamfiles_output = checkpoints.bamfile_split.get(**wildcards).output[0]
    return expand("results/{sample}/rnaseq/align/bamfiles/{file}.bam",
        sample=wildcards.sample,
        file=glob_wildcards(os.path.join(split_bamfiles_output, "{file}.bam")).file)
    

def get_rnaseq_data(wildcards):
    if rnaseq_filetype == ".bam":
        print(wildcards)
        return rnaseq_input[wildcards.sample]
    elif rnaseq_filetype == ".fastq" or rnaseq_filetype == ".fq":
        return rnaseq_input[wildcards.sample]
    else:
        print("no rnaseq data found")


# retrieve the variants
def get_variants(wildcards):
    variants = []
    if config['indel']['activate']:
        variants.append("results/wildcards.sample/rnaseq/indel/ti.indel.vcf")
        variants.append("results/wildcards.sample/rnaseq/indel/m2.indel.vcf")
    return variants



    
        


