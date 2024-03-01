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

class Counts:
    def __init__(self, countFile):
        # parse counts
        self.counts = {}
        if countFile is not None and countFile is not "":
            count_fh = open(countFile, 'r')
            lines = count_fh.readlines()
            groups = lines[0].rstrip().split('\t')[3:]
            for line in lines[1:]:
                cols = line.rstrip().split('\t')

                gene_id = cols[0]
                chrom = cols[1]

                key = (gene_id, chrom)

                if key not in self.counts:
                    self.counts[key] = {}
                    for group in groups:
                        self.counts[key][group] = float(cols[3+groups.index(group)])

            count_fh.close()

