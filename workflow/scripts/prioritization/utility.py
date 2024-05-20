from datetime import datetime

def format_output(field):
    if field == None:
        return '.'
    else:
        return field

# deter
def find_all_lowercase(seq):
    matches = []
    for i in range(len(seq)):
        if seq[i].islower():
            matches.append(i)
    return matches

# retrieve alleles from hlatyping file (e.g. mhc-I.tsv)
def get_alleles(allele_file):
    alleles = {}
    fh_alleles = open(allele_file, 'r')
    for allele in fh_alleles:
        line = allele.rstrip().split('\t')
        alleles[line[1]] = line[0]
    return alleles
    
def extract_epilens(lengths):
    """determines the lengths of the epitopes specified of the commandline 
    as these can be specific a single digits (e.g, 8,9,10) or ranges 
    (e.g., 8-10) this function extracts the single lengths"""
    epilens = []
    for lens in lengths.split(','):
        if '-' in lens:
            lower, upper = lens.split('-', 1)
            epilens.extend(range(int(lower), int(upper)+1))
        else:
            epilens.append(int(lens))
    return epilens


# determines local exon boundaries (within the transcript) - 0-based
def get_local_exon_bnds(transcript_start, exons):
    local_exons = {}
    for exon in exons.keys():
        local_exon_start = exons[exon][0] - transcript_start
        local_exon_end = exons[exon][1] - transcript_start
        local_exons[exon] = [local_exon_start, local_exon_end]
    return local_exons

# def exon_overlap(start, exons):
    # """determines if the variant overlaps with an exon"""
    # for exon in exons.keys():
        # if start >= exons[exon][0] and start <= exons[exon][0]:



def det_strand_vep(vep_strand):
    """determines the strand of the variant from the VEP annotation"""
    strand = None
    if vep_strand == '1':
        strand = "+"
    elif vep_strand == "-1":
        strand = "-"

    return strand

def transcribe(seq):
    """transcribes DNA sequence to RNA"""
    return seq.replace('T', 'U')

def revcomp(seq):
    """reverse complement a sequence"""
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([comp[base] for base in reversed(seq)])

def det_overlap(query_start, query_end, target_start, target_end):
    """determines if two regions overlap"""
    if (query_start >= target_start and query_start <= target_end or 
        query_end >= target_start and query_end <= target_end):
        return True
    else:
        return False



def get_time():
    """return the currrent time"""
    now = datetime.now()
    timestring = now.strftime("%Y-%m-%d %H:%M:%S")
    return f"[{timestring}]"




