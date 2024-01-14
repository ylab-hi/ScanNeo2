import os
import sys
import configargparse 
from pathlib import Path
import subprocess


def main():
    options = parse_arguments()
    samples = options.input.split(' ')

    counts = {}
    groups = []
    for idx, val in enumerate(options.input.split(' ')):
        groups.append(Path(val).stem.split('_counts')[0])

        fh = open(val, 'r')
        # skip header
        next(fh)
        next(fh)

        rpk_sum = 0
        for line in fh:
            val = line.rstrip().split('\t')

            gene_id = val[0]
            chrom = val[1]
            strand = val[4]

            key = (gene_id, chrom, strand)
            if key not in counts.keys():
                counts[key] = [0]*len(samples)

            count = int(val[6])
            gene_len = int(val[5])

            if count != 0:
                rpk = count/gene_len
                counts[key][idx] = count/gene_len
                rpk_sum += rpk

        fh.close()
        scale_factor = rpk_sum/1000000

        for el in counts.keys():
            counts[el][idx] = round(counts[el][idx]/scale_factor, 3)


    # generate output
    fh_output = open(options.output, 'w')
    fh_output.write('gene_id\tchrom\tstrand\t' + '\t'.join(groups) + '\n')

    for key in counts:
        fh_output.write(f'{key[0]}\t')
        fh_output.write(f'{key[1]}\t')
        fh_output.write(f'{key[2]}\t')
        fh_output.write('\t'.join(map(str, counts[key])) + '\n')

def parse_arguments():
    p = configargparse.ArgParser()

    # define command line parse_arguments
    p.add('-i', '--input', required=True, help='input file')
    p.add('-n', '--normalization', required=True, help='normalization method')
    p.add('-o', '--output', required=True, help='output file')

    # for later reference: use normal keys for DEA
    # p.add('--normal_keys', required=False, help='rna keys')

    return p.parse_args()


main()
