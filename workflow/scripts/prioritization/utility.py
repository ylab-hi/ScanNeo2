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

