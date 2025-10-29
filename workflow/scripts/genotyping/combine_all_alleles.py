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
    for infile in infiles:
        fh_in = open(infile, "r")
        for line in fh_in:
            cols = line.rstrip().split("\t")
            if len(cols) != 2:
                print(f"Invalid input file: {mhc}")
                sys.exit(1)

            source = cols[0]
            mhc = cols[1]

            # chop down the alleles to first two fields 
            # e.g., HLA-A*02:01:01 becomes HLA-A*02:01
            if len(mhc.split(':')) > 2:
                mhc = mhc.split(':')[0] + ':' + mhc.split(':')[1]


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

    if len(alleles) == 0:
        print("No valid alleles found in the input files!")
        print("Please check the input files (e.g., dnaseq, rnaseq, user-provided alleles) and try again")
        sys.exit(1)

main()
