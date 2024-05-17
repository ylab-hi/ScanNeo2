import re
import sys
import gffutils
from Bio import SeqIO

"""
This scripts modifies the header from the ENSEMBL and changes the chromosome names by adding chrs'

Usage:
    python3 modify_ensembl_header.py <fasta_file> <fasta_file_modified> <gff_file> <gff_file_modified> <nonchr>

    <nonchr> determines if non-chromosomal sequences are included as well, e.g., =0 (No) or =1 (Yes)
"""

def main():
    header_mappings = {} # header mappings
    for chr in range(1,23):
        header_mappings[str(chr)] = 'chr'+str(chr)
    header_mappings['X'] = 'chrX'
    header_mappings['Y'] = 'chrY'
    header_mappings['MT'] = 'chrMT'

    # modify fasta header
    with open(sys.argv[1], 'r') as in_fa, open(sys.argv[2], 'w') as out_fa:
        for record in SeqIO.parse(in_fa, 'fasta'):
            if record.id in header_mappings:
                record.id = header_mappings[record.id]
                SeqIO.write(record, out_fa, 'fasta')
            else:
                if sys.argv[5] == "1":
                    SeqIO.write(record, out_fa, 'fasta')

    with open(sys.argv[3], 'r') as in_gff, open(sys.argv[4], 'w') as out_gff:
        for line in in_gff:
            if line.startswith('#'):
                out_gff.write(line)
            else:
                el = line.split('\t')
                if el[0] in header_mappings:
                    el[0] = header_mappings[el[0]]
                out_gff.write('\t'.join(el))


    # # modify gtf header
    # db = gffutils.create_db(sys.argv[3], 
                            # dbfn=":memory:", 
                            # disable_infer_genes=True, 
                            # disable_infer_transcripts=True, 
                            # force=True, 
                            # keep_order=True, 
                            # merge_strategy="merge", 
                            # sort_attribute_values=True)
    # with open(sys.argv[4], 'w') as outfile:
        # for feature in db.all_features(order_by=('seqid', 'start', 'end')):
            # # Rename chromosome if it's in the mapping
            # if feature.seqid in header_mappings:
                # feature.seqid = header_mappings[feature.seqid]
            # outfile.write(str(feature) + "\n")

main()
