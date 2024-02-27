import vcfpy
import sys
from pathlib import Path

'''
This script combines multiple vcfs into one vcf. It adds a new INFO field for the group name

Usage:
python combine_vcf.py '<input_vcfs>' <method> <group> <output_vcf>
'''

def main():
    vcf = sys.argv[1]
    # extract the replicate (e.g., {rep1}_*.vcf)

    # if '_call.indel.vcf' in str(Path(vcf).name):
        # rpl = Path(vcf).name.split('_call.indel.vcf')[0]
    # elif '_somatic.short.indels.vcf' in str(Path(vcf).name):
        # rpl = Path(vcf).name.split('_somatic.short.indels.vcf')[0]
    # elif '_somatic.snvs.vcf' in str(Path(vcf).name):
        # rpl = Path(vcf).name.split('_somatic.snvs.vcf')[0]
    # elif '_fusions.vcf' in str(Path(vcf).name):
        # rpl = Path(vcf).name.split('_fusions.vcf')[0]
    # elif '_exitrons.vcf' in str(Path(vcf).name):
        # rpl = Path(vcf).name.split('_exitrons.vcf')[0]

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

    writer = vcfpy.Writer.from_path(sys.argv[4], reader.header)

    # iterate through all records
    for record in reader:
        record.INFO['SRC'] = sys.argv[2]
        record.INFO['GRP'] = sys.argv[3]
        writer.write_record(record)

main()
