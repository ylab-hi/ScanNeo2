import os
import glob
from icecream import ic
from pathlib import Path

def data_structure(data):
  config['data']['dnaseq'], filetype, readtype  = handle_seqfiles(config['data']['dnaseq'])
  config['data']['dnaseq_filetype'] = filetype
  config['data']['dnaseq_readtype'] = readtype
  config['data']['rnaseq'], filetype, readtype  = handle_seqfiles(config['data']['rnaseq'])
  config['data']['rnaseq_filetype'] = filetype
  config['data']['rnaseq_readtype'] = readtype
  return config['data']

  
def handle_seqfiles(seqdata):
  readtype = []
  filetype = []

  # iterate over replicates
  for rpl in seqdata.keys():
    files = [Path(file) for file in seqdata[rpl].split(' ')]
    ic(files)
    if len(files) == 1:  # SE
      f1_ext = get_file_extension(files[0])
      if f1_ext in ['.fq', '.fastq', '.bam']:
        seqdata[rpl] = files[0]
        filetype.append(f1_ext)
        readtype.append('SE')
      else:
        print('{} is not a valid file'.format(files[0]))
    elif len(files) == 2:  # PE
      f1_ext = get_file_extension(files[0])
      f2_ext = get_file_extension(files[1])
      # check if file extensions are the same
      if f1_ext == f2_ext:
        if(valid_paired_end(files[0], files[1])):
          seqdata[rpl] = files
          filetype.append(f1_ext)
          readtype.append('PE')
        else:
          print('files not in valid PE format')
      else:
        print('files do not have the same extension')

  # check if filetype and readtype are the same
  if all_identical(filetype) and all_identical(readtype):
    return seqdata, filetype[0], readtype[0]
  else:
    print('filetypes are not the same')
    return seqdata, None, None

# pre-processing for RNAseq data
def get_raw_reads(wildcards):
  if config['data']['rnaseq_readtype'] == 'SE':
    return config['data']['rnaseq'][wildcards.sample]
  else:  # PE
    return dict(
        zip(
            ["r1", "r2"],
            [config['data'][wildcards.seqtype][wildcards.replicate][0],
             config['data'][wildcards.seqtype][wildcards.replicate][1]],
        )
    )

# returns the reads (raw/preprocessed) for a given sample
def get_reads(wildcards):
  if config['preproc']['activate']:
    if config['data'][wildcards.seqtype+'_readtype'] == 'SE':
      return config['data'][wildcards.seqtype][wildcards.replicate]
    elif config['data'][wildcards.seqtype+'_readtype'] == 'PE':
      print("yes")
      return {"r1": "results/{sample}/{seqtype}/reads/{replicate}_preproc_r1.fq.gz",
              "r2": "results/{sample}/{seqtype}/reads/{replicate}_preproc_r2.fq.gz"}


# determines the file extension for a given file - excludes .gz
def get_file_extension(path):
  filename = path.name
  extpat = r'\.(fastq|fq|bam)(\.gz)?$'
  res = re.search(extpat, filename)
  file_ext = ''
  if res is not None:
      if res.group(0).endswith('.gz'):
          file_ext = filename[res.start():-3]
      else:
          file_ext = filename[res.start():]
  return file_ext

# check if files are a valid paired-end pair
def valid_paired_end(path1, path2):
  valid = False
  
  # only consider filename
  file1 = path1.name 
  file2 = path2.name

  # check if first file contains _R1, _1 or _fwd
  pattern = r'\_(R1|R2|1|2|fwd|rev)\.(fastq|fq){1}(\.gz)?$'
  f1_se = re.search(pattern, file1)
  f2_se = re.search(pattern, file1)

  # patterns needs to be found in both files
  if f1_se is not None and f2_se is not None:
    if file1[:f1_se.start()] == file2[:f2_se.start()]:
      valid = True
    else:
      print('{} and {} have different filestem '.format(file1, file2))
  else:
    print('{} and {} are not valid PE files'.format(file1, file2))

  return valid

# check if files in list are identical
def all_identical(l):
  if l.count(l[0]) == len(l):
    return True
  else:
    return False


config['data'] = data_structure(config['data'])
ic(config['data'])


rnaseq_filetype = ".bam"

def get_splitfastq_input(wildcards):
  if config['preproc']['activate']:
    if config['data'][wildcards.seqtype] == 'SE':
      return expand("results/{sample}/{seqtype}/preproc/reads.fq.gz", **wildcards)
    elif config['data'][wildcards.seqtype] == 'PE':  # PE
      return expand("results/{sample}/reads/{seqtype}/{replicate}_preproc_{readtype}.fq.gz",
                    readtype=["r1", "r2"],
                    seqtype=wildcards.seqtype,
                    replicate = wildcards.replicate,
                    sample = wildcards.sample)

  else:   # no pre-processing has been performed
    return rnaseq_input[wildcards.sample]


def get_splitfastq_input_PE(wildcards):
  if config['preproc']['activate']:
    return expand("results/{sample}/rnaseq/reads/inputreads_{readtype}.fq.gz",
        readtype=["r1", "r2"],
        **wildcards
    )
  else:  # no pre-processing
    return rnaseq_input[wildcards.sample]



# input for alignment w/ DNAseq data
def get_align_input_dnaseq(wildcards):
  if config['preproc']['activate']:






    if config['data']['dnaseq_readtype'] == 'SE':
      return expand("results/{sample}/dnaseq/reads/inputreads.fq.gz", **wildcards)
    else:  # PE
      return dict(
          zip(
              ["fq1", "fq2"],
              expand("results/{sample}/dnaseq/reads/inputreads_{readtype}.fq.gz",
                     readtype=["r1", "r2"],
                     **wildcards
              )
          )
      )
  else:  # no pre-processing
    if config['data']['dnaseq_readtype'] == 'SE':
      return dnaseq_input[wildcards.sample]
    else:  # PE
      return dict(
          zip(
              ["fq1", "fq2"],
              dnaseq_input[wildcards.sample]
          )
      )



def get_align_input(wildcards):
  if config['preproc']['activate']:
    if rnaseq_readtype == "se":
      return expand("results/{sample}/rnaseq/reads/inputreads.fq.gz", **wildcards)
    else:  # pe
      ic("in here")
      return dict(
          zip(
              ["fq1", "fq2"],
              expand("results/{sample}/rnaseq/reads/inputreads_{readtype}.fq.gz",
                     readtype=["r1", "r2"],
                     **wildcards
              )
          )
      )
  else:  # no pre-processing
    if rnaseq_readtype == "se":
      return rnaseq_input[wildcards.sample]
    else:  # pe
      return dict(
          zip(
              ["fq1", "fq2"],
              rnaseq_input[wildcards.sample]
          )
      )

# determine the bamfiles that contains readgroups
#def get_bams_readsgroups(wildcards):
    
#    fq1 = "results/{sample}/rnaseq/reads/{replicate}/r1/reads_{i}.fq.gz",
#    aln = "results/{sample}/rnaseq/align/{replicate}/reads_{i}.bam",


def get_readgroups_input(wildcards):
  # return only bam from STAR align
  if config['data']['rnaseq_filetype'] in ['.fq','.fastq']:
    return ["results/{sample}/rnaseq/align/{replicate}_ready.bam".format(**wildcards)]
  elif config['data']['rnaseq_readtype'] in ['.bam']:
    val = []
    val.append(config['data']['rnaseq'][wildcards.replicate][wildcards.sample])
    val.append(extend("results/{sample}/rnaseq/align/{replicate}/ready.bam",
        sample=wildcards.sample,
        replicate=wildcards.replicate))



def aggregate_alignments_fastq(wildcards):
    # make sure that all samples are processed in checkpoint - split fastq file
    checkpoint_output = checkpoints.splitfastq.get(**wildcards).output[0]
    return expand("results/{sample}/rnaseq/align/{replicate}/splt/reads_{i}.bam",
        sample=wildcards.sample,
        replicate=wildcards.replicate,
        i=glob_wildcards(os.path.join(checkpoint_output, "r1/reads_{i}.fq.gz")).i)

def aggregate_alignments_pe(wildcards):
    # make sure that all samples are processed in checkpoint - split fastq file
    checkpoint_output = checkpoints.splitfastq_pe.get(**wildcards).output[0]
    return expand("results/{sample}/rnaseq/align/bamfiles/inputreads_{i}.bam",
        sample=wildcards.sample,
        i=glob_wildcards(os.path.join(checkpoint_output, "inputreads_{i}.fq.gz")).i)

# aggregate results from STAR alignment
def aggregate_alignments(wildcards):
    # make sure that all samples are processed in checkpoint - split fastq file
    checkpoint_output = checkpoints.split_bamfile_RG.get(**wildcards).output[0]
    return expand("results/{sample}/rnaseq/align/bamfiles/{readgroup}.bam",
        sample=wildcards.sample,
        readgroup=glob_wildcards(os.path.join(checkpoint_output, "{readgroup}.bam")).readgroup)


# getting input starting from align
def get_rnaseq_data(wildcards):
    print(rnaseq_input)
    if rnaseq_filetype == ".bam":
        return rnaseq_input[wildcards.sample]
    elif rnaseq_filetype == ".fastq" or rnaseq_filetype == ".fq":
        if config['preproc']['activate']:  # preproc activated?
            if len(rnaseq_input[wildcards.sample] == 2):  # PE?
                return expand("results/{sample}/preproc/trimmed/trm.fq.gz",
                    **wildcards)
    else:
        print("no rnaseq data found")


# retrieve the variants
def get_variants(wildcards):
    variants = []
    if config['indel']['activate']:
        variants.append("results/wildcards.sample/rnaseq/indel/ti.indel.vcf")
        variants.append("results/wildcards.sample/rnaseq/indel/m2.indel.vcf")
    return variants



    
        


