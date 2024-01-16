import os
import sys

def main():
    infiles = sys.argv[1].split(' ')
    alleles = {}
    for mhc in infiles:
        fh_in = open(mhc, 'r')
        # skip header
        for line in fh_in:
            cols = line.rstrip().split('\t')
            source = cols[0]
            mhc = cols[1]

            if mhc not in alleles:
                alleles[mhc] = source
            else:
                alleles[mhc] += ',' + source
        fh_in.close()

    outfile = sys.argv[2]
    fh_out = open(outfile, 'w')
    for allele in alleles:
        fh_out.write(f'{alleles[allele]}\t{allele}\n')
    fh_out.close()


main()
