import vcfpy
import sys
from pathlib import Path

'''
This script combines multiple vcfs into one vcf. It adds a new INFO field for the group name

Usage:
python combine_vcf.py '<input_vcfs>' <method> <output_vcf>
'''

def main():
    vcfs = sys.argv[1].split(' ')
    # iterate through all vcfs
    for i,v in enumerate(vcfs):
        # extract the replicate (e.g., {rep1}_*.vcf)
        if '_sliprem.vcf' in str(Path(v).name):
            rpl = Path(v).name.split('_sliprem.vcf')[0]
        elif '_somatic.short.indels.vcf' in str(Path(v).name):
            rpl = Path(v).name.split('_somatic.short.indels.vcf')[0]
        elif '_somatic.snvs.vcf' in str(Path(v).name):
            rpl = Path(v).name.split('_somatic.snvs.vcf')[0]
        elif '_fusions.vcf' in str(Path(v).name):
            rpl= Path(v).name.split('_fusions.vcf')[0]
    

        reader = vcfpy.Reader.from_path(v)
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

        writer = vcfpy.Writer.from_path(sys.argv[3], reader.header)

        # iterate through all records
        for record in reader:
            record.INFO['GRP'] = rpl
            record.INFO['SRC'] = sys.argv[2]
            writer.write_record(record)

main()
