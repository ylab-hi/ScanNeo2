import sys
import vcfpy
import re
from pyfaidx import Fasta


def repeat_checker(s):
    if len(s) == s.count(s[0]):
        return True
    else:
        return False

def is_slippage(record):
    chrm = record.CHROM
    pos = record.POS
    alt = record.ALT
    ref = record.REF
    end = record.INFO['END']
    svlen = record.INFO['SVLEN']

    if record.INFO['SVTYPE'] == "INS":
        if repeat_checker(alt[0].value[1:]):
            pat = re.compile(rf"{alt[0].value[1] * 4}")
        else:
            pat = re.compile(rf"{alt[0].value[1:] * 4}")

    if record.INFO['SVTYPE'] == "DEL":
        if repeat_checker(ref[1:]):
            pat = re.compile(rf"{ref[1] * 4}")
        else:
            pat = re.compile(rf"{ref[1:] * 4}")
        
    match = pat.match(genome_seq[chrm][end : end + svlen * 4])
    if match:
        return True
    else:
        return False


def main():
    # read in genome sequence
    reader = vcfpy.Reader.from_path(sys.argv[2])
    writer = vcfpy.Writer.from_path(sys.argv[3], reader.header)

    for record in reader:
        if record.INFO['SVTYPE'] == 'DEL':
            # resolve ref 
            record.REF = genome_seq[record.CHROM][record.POS-1:record.INFO['END']]
            record.ALT = [vcfpy.Substitution('DEL',genome_seq[record.CHROM][record.POS-1])]

        if not is_slippage(record):
            writer.write_record(record)

genome_seq = Fasta(sys.argv[1], sequence_always_upper=True, as_raw=True)
main()











