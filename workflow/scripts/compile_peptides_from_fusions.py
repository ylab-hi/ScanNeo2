import os
import sys
import pyfaidx
import re
import configargparse
from pathlib import Path

"""
    python3 compile_peptides_from_fusions.py \
            -i '<input_fusions>' \
            -o <output_peptides> \
            -c medium \
            -p <peptide_sequences> 
            -a <annotation>
"""

def main():
    options = parse_arguments()
    # parse peptide
    pep = pyfaidx.Fasta(options.peptide)
    pep_tip = {} # list of peptide sequnece for each transcript id
    for rec in pep:
        identifier = rec.long_name
        match = re.search(r'transcript:([^.\s]+)', identifier)
        if match:
            tid = match.group(1)
            pep_tip[tid] = str(rec)

    # parse annotations
    anno = {}
    anno_fh = open(options.annotation, 'r')
    for line in anno_fh:
        if not line.startswith('#'):
            l = line.rstrip().split('\t')

            if l[2] == 'exon':
                exon_number = re.search(r'exon_number "([^.\s]+)"', l[8]).group(1)
                transcript_id = re.search(r'transcript_id "([^.\s]+)"', l[8]).group(1)

                if transcript_id not in anno:
                    anno[transcript_id] = {}
                if not exon_number in anno[transcript_id]:
                    anno[transcript_id][int(exon_number)] = [l[3],l[4]]

    
    # output file
    out_fh = open(options.output, 'w')
    write_output_header(out_fh)

    # output file for NMDs 
    out_nmd_fh = open(options.output_nmd, 'w')
    write_output_header(out_nmd_fh)


    fusions_files = options.input.split(' ')
    for fusions in fusions_files:
        # extract name of fusion file
        group = Path(fusions).stem.split('_fusions')[0]
        with open(fusions) as fh:
            next(fh)
            for line in fh:
                l = line.rstrip().split('\t')

                if l[14] == options.confidence:  # filter for selected confidence level
                    # determine if reading frame and event type are useless
                    reading_frame = l[15]
                    event_type = l[8]
                    if reading_frame == '.': # peptide sequence cannot be predicted
                        continue
                    csq = determine_consequence(reading_frame, event_type) # variant type
                    if csq is None:
                        continue
                    
                    # optional for now
                    strand = '.'
                    if l[2] != '.' and l[3] != '.':
                        strand = l[2] + '|' + l[3]
                    else:
                        continue
                   
                    # skip if peptide sequence is not available
                    peptide_seq = l[28]
                    if peptide_seq == '.':
                        continue

                    # breakpoints
                    bp1 = l[4]
                    bp2 = l[5]

                    chrom1 = bp1.split(':')[0]
                    chrom2 = bp2.split(':')[0]

                    start1 = bp1.split(':')[1]
                    start2 = bp2.split(':')[1]

                    chrom = chrom1 + '|' + chrom2
                    start = bp1.split(':')[1] + '|' + bp2.split(':')[1]
                    end = -1

                    # gene_name 
                    gene_name = l[0] + '|' + l[1]
                    gene_id = l[20] + '|' + l[21]

                    tid1 = l[22]
                    tid2 = l[23]
                    transcript_id = tid1 + '|' + tid2

                    # read (allele) depth
                    dp = f'{int(l[12])+int(l[13])}' # depth/coverage near bp1/bp2
                    mt_allele_depth = str(len(l[29]))
                    wt_allele_depth = str(int(dp) - int(mt_allele_depth))

                    # fusion transcript
                    transcript = l[27].upper()

                    # check for NMD
#                    print(f'transcript{transcript}')

                    # determine if fusion transcript is in frame 
                    transcript_breakpoint = transcript.find('|')
                    fusion_transcript = transcript.replace('|', '')

                    # print(f'breakpoint position: {break_pos}')
                    # print(f'transcript: {trscpt}')
                    # print(f'transcript2 {trscpt[break_pos:]}')

                    # consider NMD when frameshift
                    nmd = {"desc": '.', "ptc_exon_number": '.', "ptc_dist_ejc": '.', "escape_rule": '.'}
                    if 'frameshift' in csq:
                        # search for start codon in the sequence
                        start_codon = re.search(r'ATG', transcript)
                        if start_codon: # if start codon is found
                            cds = transcript[start_codon.start()+3:] # dermine coding sequence
                            cds_breakpoint = cds.find('|') # breakpoint on cds (between 1st and 2nd fusion segment)
                            cds = cds.replace('|', '') # remove breakpoint from cds

                            # search for stop_codon in cds
                            stop_pos, stop_coord = find_stop_codon(cds, cds_breakpoint, start2)

                            # check if stop codon is PTC
                            if stop_pos != -1:
                                exoninfo_tid2 = anno[tid2]
                                exons = list(exoninfo_tid2.keys())
                                exons_tid2 = list(exoninfo_tid2.keys())
                                
                                # check if stop codon is premature
                                exon_num, dist_ejc = annotate_stop_codon(exoninfo_tid2, stop_coord)
                                if exon_num != -1:
                                    nmd["ptc_exon_number"] = f'{exon_num}'
                                    if max(exons_tid2) > 1:
                                        nmd["ptc_exon_number"] += f'/{max(exons_tid2)}'
                                    nmd["ptc_dist_ejc"] = dist_ejc

                                    # check if NMD is escaped
                                    nmd_escape = check_escape(exoninfo_tid2, 
                                                              stop_coord, 
                                                              exon_num,
                                                              dist_ejc)
                                    if nmd_escape != -1:
                                        nmd["desc"] = "NMD_escaping_variant"
                                        nmd["escape_rule"] = nmd_escape


                    event = {}
                    event['chrom'] = chrom
                    event['start'] = start
                    event['end'] = end
                    event['gene_name'] = gene_name 
                    event['gene_id'] = gene_id
                    event['transcript_id'] = transcript_id
                    event['source'] = 'fusion'
                    event['group'] = group
                    event['variant_type'] = csq 

                    seg1_seq = '.'
                    if tid1 in pep_tip:
                        seg1_seq = pep_tip[tid1]
                    wt, mt, mut_pos = determine_peptide_sequence(peptide_seq, seg1_seq)
                    event['wt_subseq'] = wt
                    event['mt_subseq'] = mt
                    event['mutation_pos'] = mut_pos
                    event['vaf'] = str(-1)
                    event['wt_allele_depth'] = wt_allele_depth
                    event['mt_allele_depth'] = mt_allele_depth
                    event['read_depth'] = dp
                    event['NMD'] = nmd["desc"]
                    event['PTC_dist_ejc'] = nmd["ptc_dist_ejc"]
                    event['PTC_exon_number'] = nmd["ptc_exon_number"]
                    event['NMD_escape_rule'] = nmd["escape_rule"]

                    if event["NMD"] == "" or event["NMD"] == "NMD_escaping_variant":
                        write_output(event, out_fh)
                    else:
                        write_output(event, out_nmd_fh)

    out_fh.close()
    out_nmd_fh.close()



""" 
1. If the variant is in an intronless transcript, meaning only one exon exist in the transcript.
2. The variant falls in the first 100 coding bases in the transcript.
     vvv
  ..ES...EE..I.ES...EE.I.ES....EE.I.ES....EE 
(ES= exon_start,EE = exon_end, I = intron, v = variant location)
3. The variant location falls 50 bases upstream of the penultimate (second to the last) exon.
                           vvv
  ES...EE..I.ES...EE.I.ES....EE.I.ES....EE 
(ES= exon_start,EE = exon_end, I = intron, v = variant location)
4. The variant location  falls in the last exon of the transcript.
                                        vvvv
      ES...EE..I.ES...EE.I.ES....EE.I.ES....EE 
 (ES= exon_start,EE = exon_end, I = intron, v = variant location)

     
  
"""
def check_escape(exoninfo, stop_coord, exon_num, dist_ejc):
    exons = list(exoninfo.keys())
    last_exon = int(max(exons))

    if exon_num == 1:
        if last_exon == 1: # there is only one exon
            return 1
        elif last_exon > 1:
            exon_start = int(exoninfo[1][0])
            # check if PTC is within 
            if stop_coord - exon_start < 100:
                return 2
    if exon_num == last_exon-1:
        exon_start = int(exoninfo[last_exon-1][0])
        exon_end = int(exoninfo[last_exon-1][1])
        if exon_end - stop_coord <= 50:
            return 3
    elif exon_num == last_exon:
        return 4

    return -1


"""
    search for stop codon in cds 
    returns the index of the stop and genomic coordinate of the end of the second segment
"""
def find_stop_codon(cds, breakpoint, seg2_start):
    stop_pos = -1 # position of (detected) stop_codon --> within cds
    for i in range(0, len(cds), 3):
        codon = cds[i:i+3]
        if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
            stop_pos = i
            break

    # stop codon found before breakpoint (should be invoked)
    if stop_pos <= breakpoint:
        return -1, -1
    else:
        segment2 = stop_pos - breakpoint
        stop_coord = int(seg2_start) + segment2
        return stop_pos, stop_coord

# annotates the stop_codon (exon number and distance to exon junction)
def annotate_stop_codon(exoninfo, stop_coord):
    stop_exon = -1
    for exon_number in exoninfo:
        exon_start = int(exoninfo[exon_number][0])
        exon_end = int(exoninfo[exon_number][1])
        if stop_coord >= exon_start and stop_coord <= exon_end:
            stop_exon = exon_number
            dist_ejc = exon_end - stop_coord
            return stop_exon, dist_ejc

    # no exon information found for stop_codon
    return -1, -1


def determine_peptide_sequence(fusion_seq, wt_seq):
    # print('------')
    # print(fusion_seq)
    fusion_seg1 = fusion_seq.split('|')[0]
    fusion_seg2 = fusion_seq.split('|')[1]

    # extract sequences 
    fusion_seg1_sub = fusion_seg1.split('?')[-1]
    fusion_seg2_sub = fusion_seg2.split('?')[0].split('*')[0]
   
    # mt sequence  
    mt = (fusion_seg1_sub + fusion_seg2_sub).upper()
  
    # wt sequences
    wt = '.'
    if len(fusion_seg1_sub) > 10: # assume that the fusion sequence is long enough to contain wildtype sequence
        variants = find_all_lowercase(fusion_seg1_sub)
        if variants != []: # variants could have been found
            # extract sequence that is wildtype (in fusion - before first lowercase)
            wt_prefix = fusion_seg1_sub[:variants[0]]
        else: # no variants can be found - all of seg1 is wildtype (e.g. uppercase)
            wt_prefix = fusion_seg1_sub
        wt = extract_wildtype_sequence(wt_seq, wt_prefix, len(mt))

    return wt, mt, len(fusion_seg1_sub)+1

    
    #print(f'wildtype: {wt}, mutant: {mt}')



# extract the wildtype sequence (using ensemble info and fusion sequence)
def extract_wildtype_sequence(wt_seq, wt_prefix, mt_len):
    # search for subsequence in wildtype peptide sequences
    wt_start = wt_seq.find(wt_prefix)
#    print(f'wt_start: {wt_start}')
    if wt_start != -1: # wildtype sequence (before first variant) can be found in ensemble
        return fillup_wildtype_sequence(wt_seq[wt_start:], mt_len)
    else: # wildrtype sequence cannot be found in ensemble
        if wt_prefix.isupper():
            return fillup_wildtype_sequence(wt_prefix, mt_len)
        else: # if not all in uppercase wt cannot be determined
            return '.'

# makes sure that length of wildtype sequence equals length of mutant sequence
def fillup_wildtype_sequence(wt_seq,mt_len):
    # wildtype sequence of seg1 exceeds length of mt (take so much until length is reached)
    if len(wt_seq) >= mt_len:
        wt = wt_seq[:mt_len]
    else:
        wt = wt_seq + '$'*(mt_len-len(wt_seq))
    return wt


def find_all_lowercase(seq):
    matches = []
    for i in range(len(seq)):
        if seq[i].islower():
            matches.append(i)
    return matches


# determine consequence
def determine_consequence(reading_frame, events):
    csq = ""
    if reading_frame is not None:
        if reading_frame == "in-frame":
            csq = "inframe_"
        elif reading_frame == "out-of-frame":
            csq = "frameshift_"

        evts = events.split('/')
        for evt in evts:
            if (evt == "read-through" or
                evt == "ITD" or
                evt == "5'-5'" or
                evt == "3'-3'"):
                return None
            elif evt == "translocation":
                csq += "trs"
            elif evt == "duplication":
                csq += "dup"
            elif evt == "deletion":
                csq += "del"
            elif evt == "inversion":
                csq += "inv"


        return csq
    else:
        return None

# write
def write_output_header(stream):
    stream.write("chrom\tstart\tstop\tgene_name\tgene_id\ttranscript_id")
    stream.write("\tsource\tgroup\tvariant_type\twt_subseq\tmt_subseq")
    stream.write("\tmutation_pos\tvaf\twt_allele_depth\tmt_allele_depth")
    stream.write("\tread_depth\tNMD\tPTC_dist_ejc\tPTC_exon_number")
    stream.write("\tNMD_escape_rule")

def write_output(event, stream):
    stream.write(f'{event["chrom"]}\t')
    stream.write(f'{event["start"]}\t')
    stream.write(f'{event["end"]}\t')
    stream.write(f'{event["gene_name"]}\t')
    stream.write(f'{event["gene_id"]}\t')
    stream.write(f'{event["transcript_id"]}\t')
    stream.write(f'{event["source"]}\t')
    stream.write(f'{event["group"]}\t')
    stream.write(f'{event["variant_type"]}\t')
    stream.write(f'{event["wt_subseq"]}\t')
    stream.write(f'{event["mt_subseq"]}\t')
    stream.write(f'{event["mutation_pos"]}\t')
    stream.write(f'{event["vaf"]}\t')
    stream.write(f'{event["wt_allele_depth"]}\t')
    stream.write(f'{event["mt_allele_depth"]}\t')
    stream.write(f'{event["read_depth"]}\t')
    stream.write(f'{event["NMD"]}\t')
    stream.write(f'{event["PTC_dist_ejc"]}\t')
    stream.write(f'{event["PTC_exon_number"]}\t')
    stream.write(f'{event["NMD_escape_rule"]}\n')


def parse_arguments():
    p = configargparse.ArgParser()
    
    p.add('-i', '--input', required=True, help='Input fusion tsv file')
    p.add('-o', '--output', required=True, help='Output table')
    p.add('-c', '--confidence', required=True, choices=['high', 'medium', 'low'], help='Confidence level')  
    p.add('-l', '--log', required=False, help='Log file')
    p.add('-p', '--peptide', required=True, help='Peptide fasta file')
    p.add('-a', '--annotation', required=True, help='Annotation file')
    p.add('-n', '--output_nmd', required=True, help='Output table with NMD information')


    options = p.parse_args()
    return options

main()
    
