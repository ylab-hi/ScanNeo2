import pyfaidx
import re
import utility as ut
from intervaltree import Interval, IntervalTree

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
        print(f"{ut.get_time()} Parse reference sequences from {fastaFile}")
        self.ref = pyfaidx.Fasta(fastaFile, sequence_always_upper=True)
        print(f"{ut.get_time()} Parse transcriptome from {gtfFile}")
        self.transcriptome = self.parse_transcriptome(gtfFile) 
        print(f"{ut.get_time()} Parse exome from gtf file")
        self.exome = self.parse_exome(gtfFile)

    # 0-based coordinates and genes and transcripts
    def parse_transcriptome(self, gtfFile):
        transcriptome = {}
        fh = open(gtfFile, "r")
        gene_info = {} # store gene information
        entry = {} # entry that stores the current gene/transcript
        for line in fh:
            if not line.startswith("#"):
                el = line.rstrip().split("\t")
                chrom = el[0]
                strand = el[6]
                biotype = el[2]
                start = int(el[3])-1
                end = int(el[4])-1

                if chrom not in transcriptome:
                    transcriptome[chrom] = {}

                if biotype == "gene":
                    gene_id = self.get_attr("gene_id", el[8])
                    gene_name = self.get_attr("gene_name", el[8])
                    gene_biotype = self.get_attr("gene_biotype", el[8])


                    gene_info[gene_id] = {"gene_name": gene_name,
                                          "gene_biotype": gene_biotype,
                                          "start": start,
                                          "end": end}

                elif biotype == "transcript":
                    trx_id = self.get_attr("transcript_id", el[8])
                    trx_name = self.get_attr("transcript_name", el[8])

                    # gene id of transcript
                    trx_gid = self.get_attr("gene_id", el[8])
                    if trx_gid in gene_info: # retrieve gene information
                        trx_gene_info = gene_info[trx_gid]

                    gene_start = trx_gene_info["start"]
                    gene_end = trx_gene_info["end"]
                    transcriptome[chrom][trx_id] = {"gene_name": trx_gene_info["gene_name"],
                                                    "gene_biotype": trx_gene_info["gene_biotype"],
                                                    "gene_start": gene_start,
                                                    "gene_end": gene_end,
                                                    "transcript_strand": strand,
                                                    "transcript_name": trx_name,
                                                    "transcript_start": start,
                                                    "transcript_end": end}


        fh.close()
        return transcriptome


    def parse_exome(self, gtfFile):
        # parse exome from gtf file
        exome = {}
        current = ""
        fh = open(gtfFile, "r")
        for line in fh:
            # ignore header
            if not line.startswith("#"):
                el = line.rstrip().split("\t")
                if el[2] == "exon":
                    # extract the exon number of transcript id 
                    enum = self.get_attr("exon_number", el[8])
                    tid = self.get_attr("transcript_id", el[8])

                    start = int(el[3])-1
                    end = int(el[4])-1
                    strand = el[6]

                    if tid not in exome:
                        exome[tid] = IntervalTree()
                        current = ""

                    if current != "": # add intron
                        if strand == "+":
                            intron_start = current["end"]+1
                            intron_end = start-1

                            entry = {"strand": strand,
                                     "region": "intron",
                                     "number": f"{current['number']}-{enum}",
                                     "start": intron_start,
                                     "end": intron_end}
                        elif strand == "-":
                            intron_start = end+1
                            intron_end = current["start"]-1

                            entry = {"strand": strand,
                                     "region": "intron",
                                     "number": f"{current['number']}-{enum}",
                                     "start": intron_start,
                                     "end": intron_end}
                        exome[tid][intron_start:intron_end+1] = entry
                    
                    entry = {"strand": strand, 
                             "region": "exon",
                             "number": f"{enum}",
                             "start": start,
                             "end": end}

                    exome[tid][start:end+1] = entry
                    current = entry

        fh.close()
        return exome


    # retrieve attribute
    def get_attr(self, key, attributes):
        res = re.search(rf'{key} "([^.\s]+)"', attributes)
        if res:
            return res.group(1)
        else:
            return None



    # def parse_exome(self, gtfFile):
        # # makes more sense to have this as an interval tree

        # current = {}

        # exome = {}
        # fh = open(gtfFile, 'r')
        # for line in fh:
            # # ignore header
            # if not line.startswith('#'):
                # l = line.rstrip().split('\t')
                    # start = int(l[3])-1
                    # end = int(l[4])-1
                    # strand = l[6]

                    # if tid not in exome:
                        # exome[tid] = IntervalTree()
                        # current = ""

                    # if current != "": # add intron
                        # if strand == "+":
                            # intron_start = current["end"]+1
                            # intron_end = start-1

                            # entry = {"strand": strand,
                                     # "region": "intron",
                                     # "number": f"{current['number']}-{enum}",
                                     # "start": intron_start,
                                     # "end": intron_end}

                        # elif strand == "-":
                            # intron_start = end+1
                            # intron_end = current["start"]-1

                            # entry = {"strand": strand,
                                     # "region": "intron",
                                     # "number": f"{current['number']}-{enum}",
                                     # "start": intron_start,
                                     # "end": intron_end}

                        # # check if intron length is greater than 0
                        # if intron_start <= intron_end:
                            # exome[tid][intron_start:intron_end+1] = entry


                    # entry = {"strand": strand, 
                             # "region": "exon",
                             # "number": f"{enum}",
                             # "start": start,
                             # "end": end}

                    # exome[tid][start:end+1] = entry
                    # current = entry
                    

        # fh.close()
        # return exome

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

