import vcfpy
import sys
from pathlib import Path

'''
This script combines multiple vcfs into one vcf. It adds a new INFO field for the group name

Usage:
python combine_vcf.py '<input_vcfs>' <method> <group> <output_vcf> [<seqtype>]

seqtype records which assay produced the call. It is optional because not every
source has one: a user-supplied VCF has no assay, and leaving SEQ unset there is
more honest than inventing a value. Where a source does run on both dnaseq and
rnaseq, the same variant can be called by each, and without SEQ the resulting
rows are indistinguishable (see issue #244).
'''

def main():
    vcf = sys.argv[1]
    seqtype = sys.argv[5] if len(sys.argv) > 5 else None

    reader = vcfpy.Reader.from_path(vcf)
    reader.header.add_info_line(
        vcfpy.OrderedDict(
            [('ID', 'GRP'),
             ('Number', '1'),
             ('Type', 'String'),
             ('Description', 'Name of the group')]))
    reader.header.add_info_line(
        vcfpy.OrderedDict(
            [('ID', 'SRC'),
             ('Number', '1'),
             ('Type', 'String'),
             ('Description', 'Source of the variant')]))
    # declared even when unset, so every augmented VCF has the same header
    reader.header.add_info_line(
        vcfpy.OrderedDict(
            [('ID', 'SEQ'),
             ('Number', '1'),
             ('Type', 'String'),
             ('Description', 'Sequencing type the variant was called from')]))

    writer = vcfpy.Writer.from_path(sys.argv[4], reader.header)

    # iterate through all records
    for record in reader:
        record.INFO['SRC'] = sys.argv[2]
        record.INFO['GRP'] = sys.argv[3]
        if seqtype:
            record.INFO['SEQ'] = seqtype
        writer.write_record(record)

main()
