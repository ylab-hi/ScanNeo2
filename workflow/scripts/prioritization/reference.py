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



class Annotation:
    def __init__(self, fastaFile, gtfFile):
        self.ref = pyfaidx.Fasta(fastaFile)
        self.transcriptome = self.parse_transcriptome(gtfFile) 
        self.exome = self.parse_exome(gtfFile)

    # 0-based
    def parse_transcriptome(self, gtfFile):
        transcriptome = {}
        fh = open(gtfFile, "r")
        for line in fh:
            # ignore header
            if not line.startswith("#"):
                l = line.rstrip().split("\t")
                if l[2] == "transcript":
                    transcript_id = re.search(r'transcript_id "([^.\s]+)"', l[8]).group(1)
                    if transcript_id not in transcriptome:
                        transcriptome[transcript_id] = [int(l[3])-1, int(l[4])-1]

        fh.close()
        return transcriptome


    def parse_exome(self, gtfFile):
        exome = {}
        fh = open(gtfFile, 'r')
        for line in fh:
            # ignore header
            if not line.startswith('#'):
                l = line.rstrip().split('\t')

                if l[2] == 'exon': 
                    exon_number = re.search(r'exon_number "([^.\s]+)"', l[8]).group(1)
                    transcript_id = re.search(r'transcript_id "([^.\s]+)"', l[8]).group(1)

                    if transcript_id not in exome:
                        exome[transcript_id] = {}

                    if not exon_number in exome[transcript_id]:
                        exome[transcript_id][int(exon_number)] = [l[3],l[4]]

        fh.close()
        return exome

class Counts:
    def __init__(self, countFile):
        # parse counts
        self.counts = {}
        if countFile is not None and countFile != "":
            count_fh = open(countFile, 'r')
            lines = count_fh.readlines()
            print(lines[0].rstrip().split('\t'))
            groups = lines[0].rstrip().split('\t')[3:]
            print(groups)
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

