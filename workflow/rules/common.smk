import os
import sys
import glob
from pathlib import Path

########### CONFIG ##########
def data_structure(data):
  config['data']['dnaseq'], filetype, readtype  = handle_seqfiles(config['data']['dnaseq'])
  config['data']['dnaseq_filetype'] = filetype
  config['data']['dnaseq_readtype'] = readtype
  config['data']['rnaseq'], filetype, readtype  = handle_seqfiles(config['data']['rnaseq'])
  config['data']['rnaseq_filetype'] = filetype
  config['data']['rnaseq_readtype'] = readtype

  # abort if no data could be found
  if len(config['data']['dnaseq']) == 0 and len(config['data']['rnaseq']) == 0:
    print('No valid sequence files found')
    sys.exit(1)

  return config['data']

def handle_seqfiles(seqdata):
  readtype = []
  filetype = []

  # create new dictionary for modified information
  mod_seqdata = {}

  if seqdata is not None:
    # iterate over replicates
    for rpl in seqdata.keys():
      # make sure to ignore keys with empty values
      if seqdata[rpl] is not None:
        files = [Path(file) for file in seqdata[rpl].split(' ')]
        if len(files) == 1:  # SE
          f1_ext = get_file_extension(files[0])
          if f1_ext in ['.fq', '.fastq', '.bam']:
            mod_seqdata[rpl] = files[0]
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
              mod_seqdata[rpl] = files
              filetype.append(f1_ext)
              readtype.append('PE')
            else:
              print('files not in valid PE format')
          else:
            print('files do not have the same extension')

        # check if filetype and readtype are the same
        if all_identical(filetype) and all_identical(readtype):
          return mod_seqdata, filetype[0], readtype[0]
        else:
          print('filetypes are not the same')
          return mod_seqdata, None, None

  else:
    return mod_seqdata, None, None



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

# returns the reads (raw/preprocessed) for a given sample
def get_reads(wildcards):
  if config['preproc']['activate']:
    if config['data'][f"{wildcards.readtype}_readtype"] == 'SE':
      return config['data'][wildcards.seqtype][wildcards.replicate]
    elif config['data'][f"{wildcards.readtype}_readtype"] == 'PE':
      return {"r1": "results/{sample}/{seqtype}/reads/{replicate}_preproc_r1.fq.gz",
              "r2": "results/{sample}/{seqtype}/reads/{replicate}_preproc_r2.fq.gz"}

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

# load up the config
config['data'] = data_structure(config['data'])

########### PREPROCESSING ##########
def get_raw_reads(wildcards):
  if config['data']['rnaseq_readtype'] == 'SE':
    return dict(
        zip(
          ["sample"],
          config['data']['rnaseq'][wildcards.group]
        )
    )

# returns raw reads for a given sample
def get_qc_input(wildcards):
  return config['data'][wildcards.seqtype][wildcards.group]

def get_qc_input_fwd(wildcards):
  return config['data'][wildcards.seqtype][wildcards.group][0]

def get_qc_input_rev(wildcards):
  return config['data'][wildcards.seqtype][wildcards.group][1]

# returns the reads (raw/preprocessed) for a given sample
def get_preproc_input(wildcards):
  if config['preproc']['activate']:
    if config['data'][f"{wildcards.seqtype}_readtype"] == 'SE':
      return {
        config['data'][wildcards.seqtype][wildcards.group]
      }

    elif config['data'][f"{wildcards.seqtype}_readtype"] == 'PE':
      return {
          "sample": [config['data'][wildcards.seqtype][wildcards.group][0], 
                     config['data'][wildcards.seqtype][wildcards.group][1]]  
      }


########### HLA GENOTYPING ##########
def get_hla_flt_dna_se(wildcards):
  if config['data']['dnaseq_filetype'] == '.bam':
    return [config['data']['dnaseq'][wildcards.group]]
  else:
    if config['preproc']['activate']:
      return expand("results/{sample}/dnaseq/reads/{group}_preproc.fq.gz",
                    group = wildcards.group, sample = wildcards.sample)
    else:
      return config['data']['dnaseq'][wildcards.group]

def get_hla_flt_dna_pe(wildcards):
  if config['preproc']['activate']:
    return dict(
      zip(
        ["r1", "r2"],
        expand("results/{sample}/dnaseq/reads/{group}_{readtype}_preproc.fq.gz",
                group = wildcards.group, 
                sample = wildcards.sample, 
                readtype = ["r1","r2"])
      )
    )

  else:
    return config['data']['dnaseq'][wildcards.group]

# HLA RNA
def get_hla_flt_rna_se(wildcards):
  if config['data']['rnaseq_filetype'] == '.bam':
    return config['data']['rnaseq'][wildcards.group]
  else:
    if config['preproc']['activate']:
      return expand("results/{sample}/rnaseq/reads/{group}_preproc.fq.gz",
                    group = wildcards.group, sample = wildcards.sample)
    else:
      return config['data']['rnaseq'][wildcards.group]

def get_hla_flt_rna_pe(wildcards):
  if config['preproc']['activate']:
    return expand("results/{sample}/rnaseq/reads/{group}_{readtype}_preproc.fq.gz",
                  group = wildcards.group, 
                  sample = wildcards.sample, 
                  readtype = ["r1","r2"])
  else:
    return config['data']['rnaseq'][wildcards.group]


def hlatyping_input_DNA(wildcards):
  if config['data']['dnaseq_readtype'] == 'SE':
    return ["results/{sample}/hla/{group}_flt_dna.bam"]
  elif config['data']['dnaseq_readtype'] == 'PE':
    return {"reads": expand("results/{sample}/hla/{group}_flt_{readtype}_dna.bam",
                            sample = wildcards.sample,
                            group = wildcards.group,
                            readtype = ["r1","r2"])}

def hlatyping_input_RNA(wildcards):
  if config['data']['rnaseq_readtype'] == 'SE':
    return ["results/{sample}/hla/{group}_flt_rna.bam"]
  elif config['data']['rnaseq_readtype'] == 'PE':
    return {"reads": expand("results/{sample}/hla/{group}_flt_{readtype}_rna.bam",
                            sample = wildcards.sample,
                            group = wildcards.group,
                            readtype = ["r1","r2"])}

# returns list of hla typing results for the given sample and group
def get_alleles(wildcards):
  values = []

  if config['hlatyping']['mode'] in ['DNA', 'BOTH']:
    if config['data']['dnaseq'] is not None:
      for key in config['data']['dnaseq'].keys():
        values += expand("results/{sample}/hla/{group}_dna_result.tsv",
                           sample = wildcards.sample,
                           group = key)
    else:
      print('dnaseq data has not been specified in the config file, but specified mode for hla genotyping in config file is DNA or BOTH')


  if config['hlatyping']['mode'] in ['RNA', 'BOTH']:
    if config['data']['rnaseq'] is not None:
      for key in config['data']['rnaseq'].keys():
#        if key not in config['data']['normal']:
        values += expand("results/{sample}/hla/{group}_rna_result.tsv",
                           sample = wildcards.sample,
                           group = key)
    else:
      print('rnaseq data has not been specified in the config file, but specified mode for hla genotyping in config file is RNA or BOTH')

  if len(values) == 0:
    print('No data found. Check config file for correct specification of data and hla genotyping mode')
    sys.exit(1)

  return values



########### ALIGNMENT ##########
def get_star_input(wildcards):
  if config['data']['rnaseq_filetype'] == '.bam':
    return dict(
      zip(
        ["bam"],
        [config['data']['rnaseq'][wildcards.group]]
      )
    )

  elif config['data']['rnaseq_filetype'] == '.fq' or config['data']['rnaseq_filetype'] == '.fastq':
    if config['preproc']['activate']:
      if config['data']['rnaseq_readtype'] == 'SE':
        return expand("results/{sample}/rnaseq/preproc/reads.fq.gz", **wildcards)
      elif config['data']['rnaseq_readtype'] == 'PE':  # PE
        return dict(
            zip(
              ["fq1", "fq2"],
              expand("results/{sample}/rnaseq/reads/{group}_{readtype}_preproc.fq.gz",
                      readtype=["r1", "r2"],
                      group = wildcards.group,
                      sample = wildcards.sample)
              )
        )
      else: # no pre-processing
        return config['data']['rnaseq'][wildcards.group]

# collect the individual alignments from splitted bamfiles
def aggregate_aligned_rg(wildcards):
    # make sure that all samples are processed in checkpoint - split fastq file
    checkpoint_output = checkpoints.split_bamfile_RG.get(**wildcards).output[0]
    return expand("results/{sample}/rnaseq/align/{group}/{rg}.bam",
      sample=wildcards.sample,
      group=wildcards.group,
      rg=glob_wildcards(os.path.join(checkpoint_output, "{rg}.bam")).rg)


def get_readgroups_input(wildcards):
  # return only bam from STAR align
  if config['data'][f'{wildcards.seqtype}_filetype'] in ['.fq','.fastq']:
    return ["results/{sample}/{seqtype}/align/{group}_final_STAR.bam".format(**wildcards)]

  elif config['data'][f'{wildcards.seqtype}_filetype'] in ['.bam']:
    val = []
    
    # For DNAseq 
    if wildcards.seqtype == 'dnaseq':
      val.append(str(config['data']['dnaseq'][wildcards.group]))
    elif wildcards.seqtype == 'rnaseq':
      # needs both the raw data and star aligned bam 
      val.append(str(config['data']['rnaseq'][wildcards.group]))
      val += expand("results/{sample}/{seqtype}/align/{group}_final_STAR.bam",
          sample=wildcards.sample,
          seqtype='rnaseq',
          group=wildcards.group)
    
    return val


# input for alignment w/ DNAseq data
def get_realign_input(wildcards):
  # For DNAseq use the (raw) reads defined in config 
  if wildcards.seqtype == 'dnaseq':
    return config['data']['dnaseq'][wildcards.group]
  elif wildcards.seqtype == 'rnaseq':
    return expand("results/{sample}/{seqtype}/align/{group}_final_STAR.bam",
      sample=wildcards.sample,
      seqtype='rnaseq',
      group=wildcards.group)


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

def get_dna_align_input(wildcards):
  if config['data']['dnaseq_filetype'] in ['.fq', '.fatq']:
    if config['preproc']['activate']:
      if config['data']['dnaseq_readtype'] == 'SE':
        return expand("results/{sample}/dnaseq/reads/{group}_preproc.fq.gz", **wildcards)
      elif config['data']['dnaseq_readtype'] == 'PE':  # PE
        return expand("results/{sample}/dnaseq/reads/{group}_{readtype}_preproc.fq.gz",
                      readtype=["r1", "r2"],
                      group = wildcards.group,
                      sample = wildcards.sample)

  # is bamfile
  else:   # no pre-processing has been performed
    return config['data']['dnaseq'][wildcards.group]

########### INDEL CALLING ##########
def get_longindels(wildcards):
  indels = []
  if config['data']['dnaseq'] is not None:
    if config['indel']['mode'] in ['DNA','BOTH']:
      indels += expand("results/{sample}/{seqtype}/indel/transindel/{group}_sliprem.vcf", 
        sample=config['data']['name'], 
        seqtype='dnaseq',
        group=list(config['data']['dnaseq'].keys()))

  if config['data']['rnaseq'] is not None:
    if config['indel']['mode'] in ['RNA','BOTH']:
      indels += expand("results/{sample}/{seqtype}/indel/transindel/{group}_sliprem.vcf",
        sample=config['data']['name'], 
        seqtype='rnaseq', 
        group=list(config['data']['rnaseq'].keys()))

  return indels

def get_shortindels(wildcards):
  indels=[]
 
  if config['indel']['mode'] in ['DNA','BOTH']:
    if config['data']['dnaseq'] is not None:
      indels += expand("results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf",
        sample=config['data']['name'], 
        seqtype='dnaseq',
        group=list(config['data']['dnaseq'].keys()))
    else:
      print('dnaseq data has not been specified in the config file, but specified mode for indel calling in config file is DNA or BOTH - skipping...')
  
  if config['indel']['mode'] in ['RNA','BOTH']:
    if config['data']['rnaseq'] is not None:
      indels += expand("results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf",
        sample=config['data']['name'], 
        seqtype='rnaseq',
        group=list(config['data']['rnaseq'].keys()))
    else:
      print('rnaseq data has not been specified in the config file, but specified mode for indel calling in config file is RNA or BOTH')

  if len(indels) == 0:
    print('No data found. Check config file for correct specification of data and indel calling mode')
    sys.exit(1)

  return indels

def get_snvs(wildcards):
  indels=[]
  if config['indel']['mode'] in ['RNA','BOTH']:
    indels += expand("results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf",
      sample=config['data']['name'], 
      seqtype='rnaseq',
      group=list(config['data']['rnaseq'].keys()))
  
  if config['indel']['mode'] in ['DNA','BOTH']:
    indels += expand("results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf",
      sample=config['data']['name'], 
      seqtype='dnaseq',
      group=list(config['data']['dnaseq'].keys()))

  return indels


########### EXITRON CALLING ##########
def get_exitrons(wildcards):
  return expand("results/{sample}/rnaseq/exitron/{group}_exitron.vcf",
                sample=wildcards.sample,
                group=list(config['data']['rnaseq'].keys()))


########### GENE FUSIONS ##########
def get_fusions(wildcards):
  fusions = []
  if config['data']['rnaseq'] is not None:
    fusions += expand("results/{sample}/rnaseq/genefusion/{group}_fusions.vcf",
      sample=config['data']['name'], 
      group=list(config['data']['rnaseq'].keys()))

  return fusions


########### NEOANTIGEN PRIORIZATION ##########
def get_variants(wildcards):
    variants = []
    # indels
    if config['indel']['activate']:
      if config['indel']['type'] in ['long', 'all']:
        variants += expand("results/{sample}/annotation/long.indels.vcf",
                           sample=config['data']['name'])
      if config['indel']['type'] in ['short', 'all']:
        variants += expand("results/{sample}/annotation/somatic.short.indels.vcf",
                           sample=config['data']['name'])
        variants += expand("results/{sample}/annotation/somatic.snvs.vcf",
                           sample=config['data']['name'])

    # exitron
    if config['exitronsplicing']['activate']:
      variants += expand("results/{sample}/annotation/exitrons.vcf",
                         sample=config['data']['name'])

    # gene fusions
    if config['genefusion']['activate']:
      variants += expand("results/{sample}/annotation/fusions.vcf",
                         sample=config['data']['name'])




    return variants

