import pyfaidx
import re

class Proteome:
    def __init__(self, fastaFile):
        # parse proteome
        ref = pyfaidx.Fasta(fastaFile)
        self.proteome = {}
        for record in ref:
            identifier = record.long_name
            match = re.search(r'transcript:([^.\s]+)', identifier)
            if match:
                tid = match.group(1)
                self.proteome[tid] = str(record)


class Annotations:
    def __init__(self, gtfFile):
        # parse annotations
        anno_fh = open(gtfFile, 'r')
        self.anno = {}
        for line in anno_fh:
            # ignore header
            if not line.startswith('#'):
                l = line.rstrip().split('\t')

                if l[2] == 'exon': 
                    exon_number = re.search(r'exon_number "([^.\s]+)"', l[8]).group(1)
                    transcript_id = re.search(r'transcript_id "([^.\s]+)"', l[8]).group(1)

                    if transcript_id not in self.anno:
                        self.anno[transcript_id] = {}

                    if not exon_number in self.anno[transcript_id]:
                        self.anno[transcript_id][int(exon_number)] = [l[3],l[4]]

        anno_fh.close()

