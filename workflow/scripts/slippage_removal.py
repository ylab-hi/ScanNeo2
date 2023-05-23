import sys
import vcfpy
import re
from pyfaidx import Fasta
from icecream import ic


def repeat_checker(s):
    print(s)
    if len(s) == s.count(s[0]):
        return True
    else:
        return False

def is_slippage(record):
    chrm = record.CHROM
    pos = record.POS
    alt = record.ALT  
    ref = record.REF

    if record.INFO['SVTYPE'] == "INS":
        ic(alt[0])
        ic(alt[1])
        if repeat_checker(alt[1:]):
            pat = re.compile(rf"{alt[0]}{alt[1] * 4}")
        else:
            pat = re.compile(rf"{alt[0]}{alt[1:] * 4}")
        match = pat.match(genome_seq[chrm][pos - 1 : pos + 10])
        if match:
            return True
        else:
            return False
    elif record.INFO['SVTYPE'] == "DEL":
        ic(ref)
        if repeat_checker(ref[1:]):
            pat = re.compile(rf"{ref[0]}{ref[1] * 4}")
        else:
            pat = re.compile(rf"{ref[0]}{ref[1:] * 4}")
        match = pat.match(genome_seq[chrm][pos - 1 : pos + 10])
        if match:
            return True
        else:
            return False



def main():
    # read in genome sequence
    reader = vcfpy.Reader.from_path(sys.argv[2])



    print(genome_seq)

    for record in reader:
        if record.INFO['SVTYPE'] == 'DEL':
            # resolve ref 
            record.REF = genome_seq[record.CHROM][record.POS-1:record.INFO['END']-1]
#            print(record.REF)

        is_slippage(record)


    # writer = vcfpy.Writer.from_path(sys.argv[3], reader.header)

    # for record in reader:
        # if record.INFO["SVTYPE"] == "INS":
            # if is_slippage(record.CHROM, record.POS, record.REF, record.ALT[0].value, "INS"):
                # continue
            # else:
                # writer.write_record(record)
        # elif record.INFO["SVTYPE"] == "DEL":
            # writer.write_record(record)



#            if is_slippage(record.CHROM, record.POS, record.REF, record.ALT[0].value, "DEL"):
#                continue
#            else:
#                writer.write_record(record)
#        else:
#            writer.write_record(record)




genome_seq = Fasta(sys.argv[1], sequence_always_upper=True, as_raw=True)
main()











