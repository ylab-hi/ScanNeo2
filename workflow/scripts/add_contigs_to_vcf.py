#!/usr/bin/env python

import os
import sys
from pyfaidx import Fasta
import vcfpy
from pathlib import Path

# this a script that is required for transindel to work in the workflow 
# it adds the contigs to the vcf file

'''
    python add_contigs_to_vcf.py <input_vcf> <output_vcf> <reference_fasta>
'''

def main():
    ref = Fasta(sys.argv[3])
    vcf_reader = vcfpy.Reader.from_path(sys.argv[1])
    header = vcf_reader.header

    # add contig lines
    for chrom in ref.keys():
        header.add_contig_line(vcfpy.OrderedDict(
            [('ID', chrom),
             ('length', len(ref[chrom]))]
        ))
    
    vcf_writer = vcfpy.Writer.from_path(sys.argv[2], vcf_reader.header)
    for record in vcf_reader:
        vcf_writer.write_record(record)

main()

