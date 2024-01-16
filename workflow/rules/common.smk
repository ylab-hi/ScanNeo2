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

    for rpl in list(seqdata.keys()):
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
            #if(valid_paired_end(files[0], files[1])):
            mod_seqdata[rpl] = files
            filetype.append(f1_ext)
            readtype.append('PE')

            #else:
              #print('files not in valid PE format')
          else:
            print('files do not have the same extension')
            return mod_seqdata, None, None

        # check if filetype and readtype are the same
#        if all_identical(filetype) and all_identical(readtype):
    return mod_seqdata, filetype[0], readtype[0]

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
def get_input_hlatyping_SE(wildcards):
  if config['preproc']['activate']:
    return expand("results/{sample}/{seqtype}/reads/{group}_preproc.fq.gz",
                  sample = wildcards.sample,
                  seqtype = "dnaseq" if wildcards.type == "DNA" else "rnaseq",
                  group = wildcards.group)
  else:
    return config['data'][wildcards.type][wildcards.group]

def get_input_hlatyping_PE(wildcards):
  if config['preproc']['activate']:
    return dict(
        zip(
          ["fwd", "rev"],
          expand("results/{sample}/{seqtype}/reads/{group}_{readtype}_preproc.fq.gz",
                  sample = wildcards.sample,
                  seqtype = "dnaseq" if wildcards.type == "DNA" else "rnaseq",
                  group = wildcards.group,
                  readtype = wildcards.readtype)
          )
    )


def aggregate_mhcI_SE(wildcards):
  checkpoint_output = checkpoints.split_reads_mhcI_SE.get(**wildcards).output[0]
  return expand("results/{sample}/hla/mhc-I/genotyping/{group}_{type}_flt_SE/{no}_result.tsv",
    sample=wildcards.sample,
    group=wildcards.group,
    type=wildcards.type,
    no=glob_wildcards(os.path.join(checkpoint_output, "R_{no}.bam")).no)


def aggregate_mhcI_PE(wildcards):
  checkpoint_output = checkpoints.split_reads_mhcI_PE.get(**wildcards).output[0]
  return expand("results/{sample}/hla/mhc-I/genotyping/{group}_{type}_flt_PE/{no}_result.tsv",
    sample=wildcards.sample,
    group=wildcards.group,
    type=wildcards.type,
    no=glob_wildcards(os.path.join(checkpoint_output, "R1_{no}.bam")).no)


def get_predicted_mhcI_alleles(wildcards):
  values = []

  # routines to genotype from DNA
  if "DNA" in config['hlatyping']['MHC-I_mode']:
    if config['data']['dnaseq'] is not None:
      for key in config['data']['dnaseq'].keys():
        if key not in config['data']['normal']:
          values += expand("results/{sample}/hla/mhc-I/genotyping/{group}_{type}_{readtype}.tsv",
                           sample = wildcards.sample,
                           group = key,
                           type = "DNA",
                           readtype = config['data']['dnaseq_readtype']) # add either SE or PE
    else: # if no dnaseq data is specified, but mode is DNA or BOTH, then ignore
      print('dnaseq data has not been specified in the config file, but specified mode for hla genotyping in config file is DNA or BOTH -- will be ignored')

  # routines to genotype from RNA
  if "RNA" in config['hlatyping']['MHC-I_mode']:
    if config['data']['rnaseq'] is not None:
      for key in config['data']['rnaseq'].keys():
        if key not in config['data']['normal']:
          values += expand("results/{sample}/hla/mhc-I/genotyping/{group}_{type}_{readtype}.tsv",
                           sample = wildcards.sample,
                           group = key,
                           type = "RNA",
                           readtype = config['data']['rnaseq_readtype']) # add either SE or PE)
    else: # if no rnaseq data is specified, but mode is RNA or BOTH, then ignore
      print('rnaseq data has not been specified in the config file, but specified mode for hla genotyping in config file is RNA or BOTH -- will be ignored')


  # if alleles have been specified in the config file, add them to the list
  #if "custom" in config['hlatyping']['MHC-I_mode']:
    #if config['data']['custom']['hlatyping']['MHC-I'] is not None:
      #values.append(config['custom']['hlatyping']['MHC-I'])
  #if config['data']['custom']['hlatyping']['MHC-I'] is not None:
    #values.append(config['custom']['hlatyping']['MHC-I'])


  if len(values) == 0:
    print('No data found. Check config file for correct specification of data and hla genotyping mode')
    sys.exit(1)

  return values

def get_all_mhcI_alleles(wildcards):
  values = []

  if ("DNA" in config['hlatyping']['MHC-I_mode'] or
      "RNA" in config['hlatyping']['MHC-I_mode']):
    values += expand("results/{sample}/hla/mhc-I/genotyping/mhc-I.tsv",
                    sample = wildcards.sample)

  if "custom" in config["hlatyping"]["MHC-I_mode"]:
    values += [config["data"]["custom"]["hlatyping"]["MHC-I"]]

  if len(values) == 0:
    print('No hla data found. Check config file for correct specification of data and hla genotyping mode')
    sys.exit(1)

  return values



##### MHC CLASS I


# returns list of hla typing results for the given sample and group

###### MHC Class II #########

#def get_input_reads_classII_hlatyping(wildcards):
  #if config['preproc']['activate']:
    #return dict(
        #zip(
          #["fwd", "rev"],
          #expand("results/{sample}/{seqtype}/reads/{group}_{readtype}_preproc.fq.gz",
                 #group = wildcards.group,
                 #sample = wildcards.sample,
                 #readtype = ["r1", "r2"],
                 #seqtype = "dnaseq" if wildcards.type == "DNA" else "rnaseq")
        #)
    #)
  #else:
    #return dict(
        #zip(
          #["fwd", "rev"],
          #config['data'][wildcards.seqtype][wildcards.group]
        #)
    #)





#def get_mhcII_allels(wildcards):
  #if config['hlatyping']['MHC-II_mode'] in ['DNA', 'BOTH']:
    #if config['data']['dnaseq'] is not None:
      #if config['data']['dnaseq_readtype'] == 'PE':
        #if wildcards.
        #gT







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
                      readtype=["R1", "R2"],
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
    return expand("results/{sample}/{seqtype}/align/{group}_final_STAR.bam",
                  sample=wildcards.sample,
                  seqtype=wildcards.seqtype,
                  group=wildcards.group)

#    return ["results/{sample}/{seqtype}/align/{group}_final_STAR.bam".format(**wildcards)]

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


########### GENE EXPRESSION ##########
def get_aligned_reads(wildcards):
  val = []
  
  if wildcards.seqtype == 'dnaseq':
    val += expand("results/{sample}/{seqtype}/align/{group}_aligned_BWA.bam",
                  sample=wildcards.sample,
                  seqtype='dnaseq',
                  group=wildcards.group)
  elif wildcards.seqtype == 'rnaseq':
    val += expand("results/{sample}/{seqtype}/align/{group}_final_STAR.bam",
                  sample=wildcards.sample,
                  seqtype='rnaseq',
                  group=wildcards.group)

  return val


def get_counts(wildcards):
  counts = []

  if config['quantification']['activate']:
    if config['quantification']['mode'] in ['DNA','BOTH']:
      counts += expand("results/{sample}/{seqtype}/quantification/{group}_counts.txt",
        sample=config['data']['name'], 
        seqtype='dnaseq',
        group=list(config['data']['dnaseq'].keys()))

    if config['quantification']['mode'] in ['RNA','BOTH']:
      counts += expand("results/{sample}/{seqtype}/quantification/{group}_counts.txt",
        sample=config['data']['name'], 
        seqtype='rnaseq',
        group=list(config['data']['rnaseq'].keys()))

  return counts


########### INDEL CALLING ##########
def get_longindels(wildcards):
  indels = []
  if config['data']['dnaseq'] is not None:
    if config['indel']['mode'] in ['DNA','BOTH']:
      indels += expand("results/{sample}/{seqtype}/indel/transindel/{group}_long.indels.vcf.gz",
        sample=config['data']['name'], 
        seqtype='dnaseq',
        group=list(config['data']['dnaseq'].keys()))

  if config['data']['rnaseq'] is not None:
    if config['indel']['mode'] in ['RNA','BOTH']:
      indels += expand("results/{sample}/{seqtype}/indel/transindel/{group}_long.indels.vcf.gz",
        sample=config['data']['name'], 
        seqtype='rnaseq', 
        group=list(config['data']['rnaseq'].keys()))

  return indels

def get_shortindels(wildcards):
  indels=[]
 
  if config['indel']['mode'] in ['DNA','BOTH']:
    if config['data']['dnaseq'] is not None:
      indels += expand("results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf.gz",
        sample=config['data']['name'], 
        seqtype='dnaseq',
        group=list(config['data']['dnaseq'].keys()))
    else:
      print('dnaseq data has not been specified in the config file, but specified mode for indel calling in config file is DNA or BOTH - skipping...')
  
  if config['indel']['mode'] in ['RNA','BOTH']:
    if config['data']['rnaseq'] is not None:
      indels += expand("results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.short.indels.vcf.gz",
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
    indels += expand("results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf.gz",
      sample=config['data']['name'], 
      seqtype='rnaseq',
      group=list(config['data']['rnaseq'].keys()))
  
  if config['indel']['mode'] in ['DNA','BOTH']:
    indels += expand("results/{sample}/{seqtype}/indel/mutect2/{group}_somatic.snvs.vcf.gz",
      sample=config['data']['name'], 
      seqtype='dnaseq',
      group=list(config['data']['dnaseq'].keys()))

  return indels


########### EXITRON CALLING ##########
def get_exitrons(wildcards):
  return expand("results/{sample}/rnaseq/exitron/{group}_exitrons.vcf.gz",
                sample=wildcards.sample,
                group=list(config['data']['rnaseq'].keys()))


########### GENE FUSIONS ##########
def get_fusions(wildcards):
  fusions = []
  if config['data']['rnaseq'] is not None:
    if config['genefusion']['activate']:
      fusions += expand("results/{sample}/rnaseq/genefusion/{group}_fusions.tsv",
        sample=config['data']['name'], 
        group=list(config['data']['rnaseq'].keys()))

  return fusions

########### ALT SPLICING ##########
def get_altsplicing(wildcards):
  altsplicing = []
  if config['data']['rnaseq'] is not None:
    if config['altsplicing']['activate']:
      altsplicing += expand("results/{sample}/rnaseq/altsplicing/spladder/{group}_altsplicing.vcf.gz",
        sample=config['data']['name'], 
        group=list(config['data']['rnaseq'].keys()))
  return altsplicing


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

    return variants

def get_mhcI(wildcards):
  alleles = []
  if config['prioritization']['class'] in ['I', 'BOTH']:
    alleles += expand("results/{sample}/hla/mhc-I.tsv",
                      sample=config['data']['name'])

  return alleles

def get_mhcII(wildcards):
  alleles = []
  if config['prioritization']['class'] in ['II', 'BOTH']:
    alleles += expand("results/{sample}/hla/mhc-II.tsv",
                      sample=config['data']['name'])
  return alleles

