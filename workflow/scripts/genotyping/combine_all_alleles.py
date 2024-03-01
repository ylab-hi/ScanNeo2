import os
import sys

"""
This script combines predicted and user-provided mhc-II alleles into a
single file and also compares them with alleles from the the reference

https://help.iedb.org/hc/en-us/articles/114094151851-HLA-allele-frequencies-and-reference-sets-with-maximal-population-coverage

Usage:

python combine_all_mhcII_alleles.py '<input>' <class> <output>


class = 'mhc-I' or 'mhc-II'

"""

def parse_refset(mhc_class):
    refset = []
    with open(f"workflow/scripts/genotyping/{mhc_class}_refset.txt", "r") as fh:
        for line in fh:
            refset.append(line.rstrip())
    return refset


def main():
    infiles = sys.argv[1].split(' ')
    refset = parse_refset(sys.argv[2])

    alleles = {}
    for mhc in infiles:
        fh_in = open(mhc, "r")
        for line in fh_in:
            cols = line.rstrip().split("\t")
            source = cols[0]
            mhc = cols[1]

            if mhc in refset:
                if mhc not in alleles:
                    alleles[mhc] = source
                else:
                    alleles[mhc] += ',' + source
        fh_in.close()

    outfile = sys.argv[3]
    fh_out = open(outfile, 'w')
    for allele in alleles:
        fh_out.write(f'{alleles[allele]}\t{allele}\n')
    fh_out.close()

main()
