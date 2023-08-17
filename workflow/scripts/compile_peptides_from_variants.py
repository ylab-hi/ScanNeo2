import sys
import os
import subprocess
import re
import vcfpy
import tempfile
import configargparse
from pyfaidx import Fasta

"""
Usage:
    python variant_to_peptide.py \
            -v <input.vcf> \
            -o <outputtable> \
            -l <peptide_length> \
            -a <alleles>
"""


def main():
    options = parse_arguments()
    
    outputfile = open(options.output, 'w')
    print_header(outputfile)
    
    transcript_count = {} # counts the (different) transcripts

    input_vcfs = options.vcf.split(' ')
    for vcf in input_vcfs:
        vcf_reader = vcfpy.Reader(open(vcf, 'r'))
        csq_format = parse_csq_format(vcf_reader.header)

        transcript_count = {}
        for entry in vcf_reader:
            # FILTER (when applicable)
            if (entry.INFO['SRC'] == 'snv' or entry.INFO['SRC'] == 'short_indel'):
                if 'PASS' not in entry.FILTER:
                    continue
            
            # resolved alleles specific to vep
            alleles_vep = resolve_alleles(entry)

            chromosome = entry.CHROM
            start = entry.affected_start
            stop = entry.affected_end
            ref = entry.REF
            alts = entry.ALT

            for alt in alts:
                csq_allele = alleles_vep[str(alt.value)]
                csq_fields = parse_csq_entries(
                        entry.INFO["CSQ"], 
                        csq_format,
                        csq_allele
                )
                
                for field in csq_fields:
                    gene_name = field["SYMBOL"]
                    gene_id = field["Gene"]
                    transcript_id = field['Feature']
                    if transcript_id in transcript_count:
                        transcript_count[transcript_id] += 1
                    else:
                        transcript_count[transcript_id] = 1

                    csq = resolve_consequence(field['Consequence'])
                    if csq is None:
                        continue
                    elif csq == "frameshift":
                        if field["NMD"] != 'NMD_escaping_variant':
                            continue
                        elif field["DownstreamProtein"] == "":
                            continue

                    output_row = {
                        "chrom": entry.CHROM,
                        "start": entry.affected_start,
                        "stop": entry.affected_end,
                        "source": entry.INFO["SRC"],
                        "group": entry.INFO["GRP"],
                        "reference": entry.REF,
                        "variant": alt.value,
                        "gene_name": gene_name,
                        "transcript_id": transcript_id,
                        # "aa_change": csq_fields["Amino_acids"],
                        "gene_id": gene_id,
                        "wt_aa_seq": field["WildtypeProtein"],
                        "downstream_aa_seq": field["DownstreamProtein"],
                        "variant_type": csq,
                        "protein_position": field["Protein_position"],
                        # "seqnum": seqnum_count
                    }
                
                    if field["Amino_acids"]:
                        output_row["aa_change"] = field["Amino_acids"]
                    else:
                        output_row["aa_change"] = "NA"
                        continue
                    
                    if output_row["variant_type"] == "frameshift":
                        start_pos_var = get_variant_startpos(output_row['protein_position'])
                        wt_seq = output_row['wt_aa_seq']
                        mt_seq = wt_seq[:start_pos_var] + output_row['downstream_aa_seq']

                        # subsequence
                        wt_subseq, mt_subseq, new_start_pos_var = determine_fs_subsequences(
                            wt_seq,
                            mt_seq,
                            start_pos_var,
                        )

                        # print('###')
                        # print(f"protein position: {output_row['protein_position']}")
                        # print(f"aa change: {output_row['aa_change']}")
                        # print(f"start_pos_var: {start_pos_var}")
                        # print(f"downstream: {output_row['downstream_aa_seq']}")

                        # print(f"wt_seq: {wt_seq}")
                        # print(f"mt_seq: {mt_seq}")
                        # print(f"wt_subseq: {wt_subseq}")
                        # print(f"wt_subseq: {mt_subseq}")
                        # print(f"new protein pos: {new_start_pos_var}")

                        output_row['wt_subseq'] = wt_subseq
                        output_row['mt_subseq'] = mt_subseq
                        output_row['mutation_start'] = new_start_pos_var

                        print_row(outputfile, output_row)

                    elif (output_row['variant_type'] == 'missense' or
                          output_row['variant_type'] == 'inframe_ins' or 
                          output_row['variant_type'] == 'inframe_del'):

                        start_pos_var = get_variant_startpos(output_row['protein_position'])
                        wt_seq = output_row['wt_aa_seq']
                        wt_aa_change, mt_aa_change = determine_aa_change(output_row['aa_change'])
                        # double check that wt_aa_change is in wt_seq
                        if wt_aa_change not in wt_seq:
                            continue
                        #check for stop codons
                        wt_aa_change, wt_stop_codon = scan_stop_codon(wt_aa_change)
                        mt_aa_change, mt_stop_codon = scan_stop_codon(mt_aa_change)

                        mt_seq = wt_seq[:start_pos_var] + mt_aa_change
                        if not mt_stop_codon:
                            mt_seq += wt_seq[start_pos_var+len(wt_aa_change):]

                        # generate subsequence
                        wt_subseq, mt_subseq, new_start_pos = determine_subseq(
                            wt_seq, 
                            mt_seq,
                            start_pos_var,
                            len(wt_aa_change)
                        )


                        if 'X' in wt_subseq or 'X' in mt_subseq:
                            continue

                        if '*' in wt_subseq or '*' in mt_subseq:
                            continue

                        if len(mt_subseq) < 9:
                            continue

                        output_row['wt_subseq'] = wt_subseq
                        output_row['mt_subseq'] = mt_subseq
                        output_row['mutation_start'] = new_start_pos
                        
                        print_row(outputfile, output_row)



def determine_subseq(wt_seq, mt_seq, start_pos_var, variant_len):
    # left and right are start/end of the subsequence
    if start_pos_var < 10:
        left = 0 # start of subsequence
        new_start_pos_var = start_pos_var
    else:
        left = start_pos_var - 10
        new_start_pos_var = 10

    if start_pos_var + variant_len + 10 >= len(mt_seq):
        right = len(mt_seq)
    else:
        right = start_pos_var + variant_len + 10

    mt_subseq = mt_seq[left:right]
    if right >= len(wt_seq):
        wt_subseq = wt_seq[left:]
    else:
        wt_subseq = wt_seq[left:right]

    if wt_subseq < mt_subseq:
        for i in range(len(wt_subseq),len(mt_subseq)):
            wt_subseq += '$'



    return wt_subseq, mt_subseq, new_start_pos_var



def determine_aa_change(aa_change):
    wt_aa_change, mt_aa_change = aa_change.split('/')
    if wt_aa_change == '-':
        wt_aa_change = ''
    if mt_aa_change == '-':
        mt_aa_change = ''
    return wt_aa_change, mt_aa_change


def scan_stop_codon(sequence):
    stop_codon = False
    if '*' in sequence:
        sequence = sequence.split('*')[0]
        stop_codon = True
    if 'X' in sequence:
        sequence = sequence.split('X')[0]
        stop_codon = True

    return sequence, stop_codon



def get_variant_startpos(protein_position):
    if (protein_position is not None
        and protein_position.split('/')[0] == '-'):
        startpos = int(protein_position.split('-',1)[0])
    else:
        startpos = int(protein_position.split('-',1)[0])-1
    return startpos


def parse_csq_entries(csq_entries, csq_format, csq_allele):
    csq_format_array = csq_format.split("|")

    transcripts = []
    for entry in csq_entries:
        values = entry.split('|')
        transcript = {}
        for key, value in zip(csq_format_array, values):
            transcript[key] = value
        if transcript['Allele'] == csq_allele:
            transcripts.append(transcript)

        return transcripts


def parse_csq_format(vcf_header):
    if vcf_header.get_info_field_info('CSQ') is None:
        sys.exit("Failed to extract format string form info description for tag (CSQ)")
    else:
        csq_header = vcf_header.get_info_field_info('CSQ')
        format_pattern = re.compile("Format: (.*)")
        match = format_pattern.search(csq_header.description)
        return match.group(1)

def determine_fs_subsequences(wt_seq, mt_seq, start_pos_var):
    subseqs = {}
    if start_pos_var < 10:
        subseq_start = 0
    else:
        subseq_start = start_pos_var - 10

    wt_subseq = wt_seq[subseq_start:]
    mt_subseq = mt_seq[subseq_start:]

    # make sure wt and mt subsequence of of same length
    if len(wt_subseq) > len(mt_subseq):
        wt_subseq = wt_subseq[:len(mt_subseq)]
    else:
        for i in range(len(wt_subseq),len(mt_subseq)):
            wt_subseq += '$'

    # determine (new) variant start within subsequence
    new_start_pos_var = start_pos_var - subseq_start

    return wt_subseq, mt_subseq, new_start_pos_var


def resolve_alleles(entry):
    alleles = {}
    for alt in entry.ALT:
        alt = str(alt.value)
        # most likely an SNV
        if alt[0:1] != entry.REF[0:1]:
            alleles[alt] = alt
        elif alt[1:] == "":
            alleles[alt] = "-"
        else:
            alleles[alt] = alt[1:]

    return alleles

def resolve_consequence(consequence_string):
    consequences = {
        consequence.lower() for consequence in consequence_string.split("&")
    }
    if "start_lost" in consequences:
        consequence = None
    elif "frameshift_variant" in consequences:
        consequence = "frameshift"
    elif "missense_variant" in consequences:
        consequence = "missense"
    elif "inframe_insertion" in consequences:
        consequence = "inframe_ins"
    elif "inframe_deletion" in consequences:
        consequence = "inframe_del"
    else:
        consequence = None
    
    return consequence

def print_header(outputfile):
    # keys for correct order
    keys = ["chrom","start", "stop","gene_name",
            "gene_id","transcript_id","source",
            "group","variant_type", "wt_subseq",
            "mt_subseq", "mutation_pos"
    ]
    
    # print header
    for key in keys:
        outputfile.write(key + '\t')
    outputfile.write('\n')


def print_row(outputfile, row):
    outputfile.write(f"{row['chrom']}\t{row['start']}\t{row['stop']}")
    outputfile.write(f"\t{row['reference']}\t{row['variant']}")
    outputfile.write(f"\t{row['gene_name']}\t{row['gene_id']}")
    outputfile.write(f"\t{row['transcript_id']}\t{row['source']}")
    outputfile.write(f"\t{row['group']}\t{row['variant_type']}")
    outputfile.write(f"\t{row['wt_subseq']}\t{row['mt_subseq']}")
    outputfile.write(f"\t{row['mutation_start']}\n")



def parse_arguments():
    p = configargparse.ArgParser()
    
    p.add('-v', '--vcf', required=True, help='Input vcf files')
    p.add('-o', '--output', required=True, help='Output table')
    p.add('-l', '--log', required=False, help='Log file')

    options = p.parse_args()

    return options

main()
