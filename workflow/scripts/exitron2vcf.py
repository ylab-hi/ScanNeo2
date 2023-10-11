#!/usr/bin/env python

import os
import sys
import argparse
from pyfaidx import Fasta
import vcfpy
from pathlib import Path

'''
    python exitron2vcf2.py <input_exitron> <output_vcf> <reference_fasta>
'''

def add_info_field(header, key, number, dtype, desc):
    info = vcfpy.OrderedDict([
        ('ID', key),
        ('Number', number),
        ('Type', dtype),
        ('Description', desc)
    ])
    header.add_info_line(info)
    return header

def create_header(samplename, ref):
    samples = vcfpy.SamplesInfos([samplename])
    # extract samples from filename
    header = vcfpy.Header(
            lines = [
                    vcfpy.HeaderLine(key="fileformat", value="VCFv4.2"),
            ],
            samples=samples)

    # add fields to header
    header = add_info_field(header, 'DP', 1, 'Integer', 'Total depth of junction')
    header = add_info_field(header, 'SplicedSite', '1', 'String', 'Spliced site')
    header = add_info_field(header, 'STRAND', '1', 'String', 'Junction Strand')
    header = add_info_field(header, 'END', '1', 'Integer', 'End position of the structural variant')
    header = add_info_field(header, 'AO', '1' , 'Integer', 'Number of reads supporting the exitron')
    header = add_info_field(header, 'SVTYPE', '1', 'String', 'Type of structural variant')
    header = add_info_field(header, 'SVLEN', '1', 'Integer', 'Length of the structural variant')
    header = add_info_field(header, 'PSO', '1', 'Float', 'Percentage spliced-out')
    header = add_info_field(header, 'PSI', '1', 'Float', 'Percentage spliced-in')
    header = add_info_field(header, 'GeneName', '1', 'String', 'Gene name')
    header = add_info_field(header, 'GeneID', '1', 'String', 'Gene ID')
    header = add_info_field(header, 'MAPQ', '1', 'Integer', 'Median mapping quality of paired-ends')
    header = add_info_field(header, 'SR', '1', 'Integer', 'Number of split reads supporting the exitron')
    header = add_info_field(header, 'SRQ', '1', 'Float', 'Split-read consensus alingment quality')
    header = add_info_field(header, 'CONSENSUS', '1', 'String', 'Split-read consensus sequence')
    header = add_info_field(header, 'CE', '1', 'Float', 'Consensus sequence entropy')
    header = add_info_field(header, 'CT', '1', 'String', 'Paired-end signature induced connection type')
    header = add_info_field(header, 'IMPRECISE', '0', 'Flag', 'Imprecise structural variation')
    header = add_info_field(header, 'PRECISE', '0', 'Flag', 'Precise structural variation')
    header = add_info_field(header, 'SVMETHOD', '1', 'String', 'Type of approach used to detect SV')
    header = add_info_field(header, 'INSLEN', '1', 'Integer', 'Predicted length of the insertion')
    header = add_info_field(header, 'HOMLEN', '1', 'Integer', 'Predicted microhomology length using a max. edit distance of 2')
    header.add_format_line(vcfpy.OrderedDict(
        [('ID','GT'),
         ('Number',1),
         ('Type', 'String'),
         ('Description', 'Genotype')
         ]))

    # add contig lines
    for chrom in ref.keys():
        header.add_contig_line(vcfpy.OrderedDict(
            [('ID', chrom),
             ('length', len(ref[chrom]))]
        ))


    return header


def main():
    seq = Fasta(sys.argv[3])

    samplename = Path(sys.argv[1]).stem

    header = create_header(samplename,seq)
    writer = vcfpy.Writer.from_path(sys.argv[2], header)

    with open(sys.argv[1]) as fh:
        next(fh)
        for line in fh:
            el = line.rstrip().split('\t')
            chrom = el[0]
            pos = int(el[1])
            end = int(el[2])
            _alt = str(seq[chrom][pos-1:pos])
            _ref = str(seq[chrom][pos-1:end-1])
            idx = el[3]
            ao = int(el[4])
            strand = el[5]
            gene_name = el[6]
            length = int(el[7])
            spliced_site = el[8]
            gene_id = el[9].split('.')[0]
            pso = float(el[10])
            psi = float(el[11])
            dp = int(el[12])

            rec = vcfpy.Record(
                    CHROM=chrom,
                    POS=pos,
                    ID=[idx],
                    REF=_ref,
                    ALT=[vcfpy.Substitution(type_='DEL', value=_alt)],
                    QUAL='.',
                    FILTER=['PASS'],
                    INFO = vcfpy.OrderedDict(
                        [('SVTYPE', 'DEL'),
                         ('END', end-1),
                         ('DP', dp),
                         ('STRAND', strand),
                         ('SplicedSite', spliced_site),
                         ('GeneName', gene_name),
                         ('SVLEN', int('-'+str(length))),
                         ('PSO', pso),
                         ('PSI', psi)]
                    ),
                    FORMAT=['GT'],
                    calls = [vcfpy.Call(sample=samplename, data={'GT': '0/1'})]
            )
            writer.write_record(rec)

main()    
