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
                    
    # retrieve the different epitope_lens
    epi_lo, epi_up = options.peptide_length.split('-')

    outputfile = open(options.output, 'w')
    wt_fasta = open('workflow/scripts/wt.fasta','w')
    mt_fasta = open('workflow/scripts/mt.fasta','w')

    # keys for correct order
    keys = ["chrom","start", "stop","gene_name",
            "gene_id","transcript_id","source",
            "group","variant_type","wt_epitope_seq",
            "wt_epitope_ic50","wt_epitope_rank",
            "mt_epitope_seq", "mt_epitope_ic50",
            "mt_epitope_rank", "allele", "mutation_pos"
    ]

    # print header
    for key in keys:
        outputfile.write(key + '\t')
    outputfile.write('\n')

    output_entries = []

    # parse alleles
    alleles_dict = {}
    alleles_fh = open(options.alleles, 'r')
    for al in alleles_fh:
        group, alleles = al.split('\t',1)
        alleles_dict[group] = []
        for x in alleles.split('\t'):
            alleles_dict[group].append(x.rstrip())

    transcript_count = {} # counts the (different) transcripts

    all_subseqs = {}
    for epitope_length in range(int(epi_lo), int(epi_up)+1):
        all_subseqs[epitope_length] = {}
    seqnum_count = 1

    # for epitope_length in range(int(epi_lo), int(epi_up)+1):
        # wt_fasta = open('workflow/scripts/wt.fa','w')
        # mt_fasta = open('workflow/scripts/mt.fa','w')

    wt_fasta = open('workflow/scripts/wt.fa', 'w')
    mt_fasta = open('workflow/scripts/mt.fa', 'w')

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

                if len(csq_fields) == 0:
                    continue
                else:
                    csq_fields = csq_fields[0]

                gene_name = csq_fields["SYMBOL"]
                gene_id = csq_fields["Gene"]
                transcript_id = csq_fields['Feature']
                if transcript_id in transcript_count:
                    transcript_count[transcript_id] += 1
                else:
                    transcript_count[transcript_id] = 1

                csq = resolve_consequence(csq_fields['Consequence'])
                if csq is None:
                    continue
                elif csq == "frameshift":
                    if csq_fields["NMD"] != 'NMD_escaping_variant':
                        continue
                    elif csq_fields["DownstreamProtein"] == "":
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
                    "aa_change": csq_fields["Amino_acids"],
                    "gene_id": gene_id,
                    "wt_aa_seq": csq_fields["WildtypeProtein"],
                    "downstream_aa_seq": csq_fields["DownstreamProtein"],
                    "variant_type": csq,
                    "protein_position": csq_fields["Protein_position"],
                    "seqnum": seqnum_count
                }

                if csq_fields["Amino_acids"]:
                    output_row["aa_change"] = csq_fields["Amino_acids"]
                else:
                    output_row["aa_change"] = "NA"

                if output_row["variant_type"] == "frameshift":
                    
                    # position at which variations occurs
                    start_pos_var = get_variant_startpos(output_row['protein_position'])
                    wt_seq = output_row['wt_aa_seq']
                    mt_seq = wt_seq[:start_pos_var] + output_row['downstream_aa_seq']

                    # subsequence
                    wt_subseq, mt_subseq, new_start_pos_var = determine_fs_subsequences(
                        wt_seq,
                        mt_seq,
                        epitope_length,
                        start_pos_var,
                    )

                    output_row['wt_subseq'] = wt_subseq
                    output_row['mt_subseq'] = mt_subseq
                    output_row['new_var_start_pos'] = new_start_pos_var

                    # remove dollar signs & and lengths at least to worj with mhci
                    wt_subseq_stripped = wt_subseq.split('$',1)[0]
                    if len(wt_subseq_stripped) >= 8:
                        wt_fasta.write(f">{seqnum_count}\n{wt_subseq.split('$',1)[0]}\n")
                    if len(mt_subseq) >= 8:
                        mt_fasta.write(f">{seqnum_count}\n{mt_subseq}\n")

                    output_entries.append(output_row)
                    
                    seqnum_count += 1

                elif (output_row['variant_type'] == "missense" or
                      output_row['variant_type'] == 'inframe_ins' or 
                      output_row['variant_type'] == 'inframe_del'):
                    # print(output_row)
                    start_pos_var = get_variant_startpos(output_row['protein_position'])
                    wt_seq = remove_stop_codons(output_row['wt_aa_seq'])
                    mt_seq = remove_stop_codons(wt_seq[:start_pos_var] + output_row['variant'] + wt_seq[start_pos_var+len(output_row['variant']):])
                    
                    # generate subsequence
                    wt_subseq, mt_subseq, new_start_pos = determine_subseq(
                        wt_seq,
                        mt_seq,
                        start_pos_var,
                        len(output_row['variant'])
                    )

                    # print(f"wt {wt_seq}")
                    # print(f"mt {mt_seq}")

                    output_entries.append(output_row)

                    seqnum_count += 1
                    

    # no significant entries found.. abort
    if len(output_entries) == 0:
        sys.exit()

    wt_fasta.close()
    mt_fasta.close()
    wt_affinity = calc_peptide_binding(alleles_dict,
                         'workflow/scripts/wt.fa', 
                         list(range(int(epi_lo),int(epi_up)+1)),'wt')

    mt_affinity = calc_peptide_binding(alleles_dict,
                         'workflow/scripts/mt.fa', 
                         list(range(int(epi_lo),int(epi_up)+1)),'mt')

    print(mt_affinity)


    immuno_seqs = open('workflow/scripts/im.txt','w')
    for entry in output_entries:
        final_result = entry
        seqnum = final_result['seqnum']

        wt = None
        # search for epitopes with high binding affinity
        if seqnum in wt_affinity:
            wt = wt_affinity[seqnum]
        else:
            final_result['wt_epitope_ic50'] = 'NA'
            final_result['wt_epitope_rank'] = 'NA'
        
        if seqnum in mt_affinity:
            mt = mt_affinity[seqnum]
        else:
            continue

        for epitope in mt.keys():
            # determine by mhc_i / convert to 0-based
            start_pos_in_subseq = mt[epitope][1]-1
            end_pos_in_subseq = mt[epitope][2]-1
    
            # check if mutation is part of the subsequence (within or upstream)
            if final_result['new_var_start_pos'] >= end_pos_in_subseq:
                continue
            elif final_result['new_var_start_pos'] <= start_pos_in_subseq:
                final_result['mutation_position'] = 0
            else:
                final_result['mutation_position'] = final_result['new_var_start_pos']

            final_result['mt_epitope_seq'] = epitope
            final_result['allele'] = mt[epitope][0]
            final_result['mt_epitope_ic50'] = mt[epitope][3]
            final_result['mt_epitope_rank'] = mt[epitope][4]

            # search for corresponding WT
            startpos_epitope_subseq = final_result['mt_epitope_seq'].find(epitope)
            final_result['wt_epitope_seq'] = final_result['wt_subseq'][startpos_epitope_subseq:len(epitope)]

            if wt is not None:
                if final_result['wt_epitope_seq'] in wt.keys():
                    final_result['wt_epitope_ic50'] = wt[final_result['wt_epitope_seq']][3]
                    final_result['wt_epitope_rank'] = wt[final_result['wt_epitope_seq']][4]
                else:
                    final_result['wt_epitope_ic50'] = 'NA'
                    final_result['wt_epitope_rank'] = 'NA'

            print(final_result)
            print_entries(final_result,outputfile)
            # write back to txt for immunogenicity
            immuno_seqs.write(f"{final_result['mt_epitope_seq']}\n")

           # print(final_result)
    
    outputfile.close()  # close initial output file
    immuno_seqs.close()
    immunogenecity_values = ['immunogenicity'] + calc_immunogenecity('workflow/scripts/im.txt')
    # add immunogenicity scores to result table
    with open(options.output, 'r+') as file:
        lines = file.readlines()
        file.seek(0,0)
        for i,v in enumerate(lines):
            row = v.rstrip().split('\t') + [str(immunogenecity_values[i])]
            file.write('\t'.join(row)+'\n')



def print_entries(row, out):
    out.write(f"{row['chrom']}\t")
    out.write(f"{row['start']}\t")
    out.write(f"{row['stop']}\t")
    out.write(f"{row['gene_name']}\t")
    out.write(f"{row['gene_id']}\t")
    out.write(f"{row['transcript_id']}\t")
    out.write(f"{row['source']}\t")
    out.write(f"{row['group']}\t")
    out.write(f"{row['variant_type']}\t")
    out.write(f"{row['wt_epitope_seq']}\t")
    out.write(f"{row['mt_epitope_seq']}\t")
    out.write(f"{row['wt_epitope_ic50']}\t")
    out.write(f"{row['mt_epitope_ic50']}\t")
    out.write(f"{row['allele']}\t")
    out.write(f"{row['mutation_position']}\t")
    out.write(f"\n")



def determine_subseq(wt_seq, mt_seq, var_start_pos, var_len):
    if var_start_pos < 10:
        adj_start_pos = 0
    else:
        adj_start_pos = var_start_pos - 10

    if len(mt_seq)-(var_start_pos + 1) > var_len+10:
        mt_subseq = mt_seq[adj_start_pos:var_start_pos+var_len+10]
        # what if wt is shorter: probably insertion
        if len(wt_seq[adj_start_pos:]) < len(mt_subseq):
            wt_subseq = wt_seq[adj_start_pos:]
            for i in range(len(wt_subseq),len(mt_subseq)):
                wt_subseq += '$'
        else:
            wt_subseq = wt_seq[adj_start_pos:var_start_pos+var_len+10]
    else:
        mt_subseq = mt_seq[adj_start_pos:]
        wt_subseq = wt_seq[adj_start_pos:len(mt_subseq)]

    new_start_pos_var = var_start_pos = adj_start_pos

    return wt_subseq, mt_subseq, new_start_pos_var



def determine_fs_subsequences(wt_seq, mt_seq, epitope_length, start_pos_var):
    subseqs = {}
    if epitope_length > start_pos_var:
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


        # # subseq_start = start - (epitope_length - 1)
    # # subseq_end = start + len(wt_seq) - 1 

    # mt_subseq = mt_seq[subseq_start:]

    # if wt_seq >= mt_seq:
        # wt_subseq = wt_seq[subseq_start:subseq_end]
    # else:
        # wt_subseq = wt_seq[subseq_start:]

    # return wt_subseq, mt_subseq


    # print(wildtype_seq)
    # print(mutant_seq)
    
    # print(start)
    # print(end)

    # print(f"subseq_start: {subseq_start}")
    # print(f"subseq_end: {subseq_end}")

    # for i in range(subseq_start, len(mutant_seq)-epitope_length+1):
        # if ((i >= start and i <= end) or 
        # (i+epitope_length-1 >= start and i+epitope_length-1 <= end)):
            # # print(f"i: {i}")
            # # print(f"i: {i+epitope_length}")

            # # print(f"wt: {wildtype_seq[i:i+epitope_length]}")
            # # print(f"mt: {mutant_seq[i:i+epitope_length]}")

            # wt_subseq = wildtype_seq[i:i+epitope_length]
            # mt_subseq = mutant_seq[i:i+epitope_length]

            # # only consider peptide if not already detected
            # if mt_subseq not in subseqs:
                # subseqs[mt_subseq] = (wt_subseq, mt_subseq)

    # return subseqs


# calculate binding affinity
def calc_peptide_binding(alleles, fa_file, epilens, wt_mt ):
    binding_affinity = {}
    for algroup in alleles.keys():
        allele_list = alleles[algroup]
        for al in allele_list:
            al_mod = 'HLA-'+al
            if valid_alleles(al_mod):
                for epilen in epilens:
                    call = ['python', 
                        'workflow/scripts/mhc_i/src/predict_binding.py', 
                        'netmhcpan', 
                        al_mod, 
                        str(epilen),
                        fa_file]
                    result = subprocess.run(call,
                        stdout = subprocess.PIPE,
                        universal_newlines = True 
                    )
                    predictions = result.stdout
                    line_by_line = predictions.rstrip().split('\n')
                    for line in line_by_line[1:]:
                        entries = line.split('\t')
                        if wt_mt == 'mt':
                            if float(entries[8]) >= 500:
                                continue
                        
                        if int(entries[1]) not in binding_affinity:
                            binding_affinity[int(entries[1])] = {}

                        if entries[5] not in binding_affinity[int(entries[1])]:
                            # allele, start, end, ic50, rank
                            binding_affinity[int(entries[1])][entries[5]] = (entries[0], int(entries[2]), int(entries[3]), float(entries[8]), float(entries[9]))


    return binding_affinity



def calc_immunogenecity(seq):
    print(seq)
    result = subprocess.run(
        ['python', 
         'workflow/scripts/immunogenicity/predict_immunogenicity.py', 
         str(seq)],
        stdout = subprocess.PIPE,
        universal_newlines = True 
    )

    res = result.stdout.rstrip().split('\n')[4:]
    scores = [float(item.split(',')[2]) for item in res if len(item.split(',')) >= 3]
    return scores



#
def determine_variant_length(reference, variant):
    # variation
    if len(reference) == len(variant):
        return 1
    else:
        # deletion: returns negative value
        # insertion: returnr positive value
        return len(variant) - len(reference)





    


def get_variant_startpos(protein_position):
    if (protein_position is not None
        and protein_position.split('/')[0] == '-'):
        startpos = int(protein_position.split('-',1)[0])
    else:
        startpos = int(protein_position.split('-',1)[0])-1
    return startpos


def remove_stop_codons(seq):
    if '*' in seq:
        seq = seq.split('*')[0]
    elif 'X' in seq:
        seq = seq.split('X')[0]
    return seq



def get_frameshift_subseq(pos, wt_seq, pep_seq_len, output_row):
    flanking_seq_len = determine_flanking_seq_len(
            len(wt_seq), pep_seq_len, output_row
    )
    print(f"flanking_seq_len: {flanking_seq_len}")
    if pos < flanking_seq_len:
        start_pos = 0
    else:
        start_pos = int(pos - flanking_seq_len)
    wt_subseq_stop = int(pos + flanking_seq_len)
    mt_subseq_stop = int(pos)

    wt_subseq = wt_seq[start_pos:wt_subseq_stop]
    mt_subseq = wt_seq[start_pos:wt_subseq_stop]

    return (start_pos, wt_subseq, mt_subseq)


def get_wt_subseq(pos, wt_seq, wt_seq_len, pep_seq_len, output_row):
    flanking_seq_len = int(
            determine_flanking_seq_len(
                len(wt_seq), pep_seq_len, output_row
            )
    )
    pep_seq_len = min(
            2 * flanking_seq_len + wt_seq_len,
            len(wt_seq)
    )

    if distance_from_start(pos, wt_seq) < flanking_seq_len:
        wt_subseq = wt_seq[:pep_seq_len]
        mt_pos = pos 
    elif distance_from_end(pos, wt_seq) < flanking_seq_len:
        start_pos = len(wt_seq) - pep_seq_len
        wt_subseq = wt_seq[start_pos:]
        mt_pos = pep_seq_len - distance_from_end(pos, wt_seq) - 1
    elif (distance_from_start(pos, wt_seq) >= flanking_seq_len and
          distance_from_end(pos, wt_seq) >= flanking_seq_len):
        start_pos = pos - flanking_seq_len
        end_pos = start_pos + pep_seq_len
        wt_subseq = wt_seq[start_pos:end_pos]
        mt_pos = flanking_seq_len
    else:
        sys.exit('ERROR reading wildtype subsequence')

    return mt_pos, wt_subseq


def distance_from_start(pos, seq):
    return pos

def distance_from_end(pos, seq):
    return len(seq) - 1 - pos


def determine_flanking_seq_len(wt_seq_len, pep_seq_len, output_row):
    # actual possible peptide length
    corr_pep_seq_len = determine_pep_seq_len(wt_seq_len, pep_seq_len, output_row)
    print(f"possible peptide length: {corr_pep_seq_len}")
    if corr_pep_seq_len % 2 == 0:
        return (corr_pep_seq_len - 2) / 2
    else:
        return (corr_pep_seq_len - 1) / 2


def determine_pep_seq_len(wt_seq_len, pep_seq_len, output_row):
    corr_pep_seq_len = pep_seq_len
    if wt_seq_len < pep_seq_len:
        corr_pep_seq_len = wt_seq_len
        print("WT len shorter than desired peptide len")
    return corr_pep_seq_len



# this function converts the alleles in an vcfentry to vep specific ones
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


# determines the downstream effect of the variant
def determine_mt_aa_seq(csq, wt_seq, pos, aa_change):
#    if csq == 'frameshift_variant':

    aa_wt_change = aa_change.split('/')[0]
    aa_wt = '' if aa_wt_change == '-' else aa_wt_change
    aa_mt_change = aa_change.split('/')[1]
    aa_mt = '' if aa_mt_change == '-' else aa_mt_change

    aa_pos = pos.split('-')
    if len(aa_pos) == 1:  # single amino acid change
        start_pos = int(aa_pos[0])
        end_pos = start_pos+1
        mt_seq = wt_seq[0:start_pos-1] + aa_mt + wt_seq[end_pos:]
        end_new = start_pos + len(aa_mt)
    else:
        start_pos = int(aa_pos[0])
        end_pos = int(aa_pos[1])
        mt_seq = wt_seq[0:start_pos-1] + aa_mt + wt_seq[end_pos:]
        end_new = start_pos + len(aa_mt)-1

    return mt_seq


def iedb_call(allele, peptide, size):
    with tempfile.TemporaryDirectory() as tmpdirname:
        print('created temporary directory', tmpdirname)

        fasta = tempfile.NamedTemporaryFile()

        print(fasta)



def determine_subpeptide(mt_seq, pos, aa_change):
    # check if the variant is a single amino acid change
    start = int(pos.split('-')[0])  # start position of variants
    wt_aa = aa_change.split('/')[0]  # wildtype amino acid
    mt_aa = aa_change.split('/')[1]  # mutant amino acid

    if mt_aa == '-':
        end = start
    else:
        end = start + len(mt_aa)

    if wt_aa == '-':
        end = start + len(mt_aa) - 1
    else:
        end = start + len(mt_aa)

    if start - 10 <= 0:
        left = 1
    else:
        left = start - 10

    if end + 10 >= len(mt_seq):
        right = len(mt_seq)
    else:
        right = end + 10

    return mt_seq[left-1:right-1]


def simplify_indel_allele(ref, alt):
    while len(ref) > 0 and len(alt.value) > 0 and ref[-1] == alt.value[-1]:
        ref = ref[0:-1]
        alt_new = alt.value[0:-1]
    while len(ref) > 0 and len(alt.value) > 0 and ref[0] == alt.value[0]:
        ref = ref[1:]
        alt_new = alt.value[1:]
    return ref, alt_new

            


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


def parse_csq_entries_for_allele(csq_entries, csq_format, csq_allele):
    csq_format_array = csq_format.split("|")

    transcripts = []
    for entry in csq_entries:
        values = entry.split("|")
        transcript = {}
        for key, value in zip(csq_format_array, values):
            transcript[key] = value
        if transcript["Allele"] == csq_allele:
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


def parse_arguments():
    p = configargparse.ArgParser()
    
    p.add('-v', '--vcf', required=True, help='Input vcf files')
    p.add('-o', '--output', required=True, help='Output table')
    p.add('-p', '--peptide_length', required=True, help='Peptide length')
    p.add('-a', '--alleles', required=True, help='Alleles')
    p.add('-l', '--log', required=False, help='Log file')

    options = p.parse_args()

    return options



def valid_alleles(allele):
    allele_list = [
        "HLA-A*01:01",
                "HLA-A*01:02",
                "HLA-A*01:03",
                "HLA-A*01:06",
                "HLA-A*01:07",
                "HLA-A*01:08",
                "HLA-A*01:09",
                "HLA-A*01:10",
                "HLA-A*01:12",
                "HLA-A*01:13",
                "HLA-A*01:14",
                "HLA-A*01:17",
                "HLA-A*01:19",
                "HLA-A*01:20",
                "HLA-A*01:21",
                "HLA-A*01:23",
                "HLA-A*01:24",
                "HLA-A*01:25",
                "HLA-A*01:26",
                "HLA-A*01:28",
                "HLA-A*01:29",
                "HLA-A*01:30",
                "HLA-A*01:32",
                "HLA-A*01:33",
                "HLA-A*01:35",
                "HLA-A*01:36",
                "HLA-A*01:37",
                "HLA-A*01:38",
                "HLA-A*01:39",
                "HLA-A*01:40",
                "HLA-A*01:41",
                "HLA-A*01:42",
                "HLA-A*01:43",
                "HLA-A*01:44",
                "HLA-A*01:45",
                "HLA-A*01:46",
                "HLA-A*01:47",
                "HLA-A*01:48",
                "HLA-A*01:49",
                "HLA-A*01:50",
                "HLA-A*01:51",
                "HLA-A*01:54",
                "HLA-A*01:55",
                "HLA-A*01:58",
                "HLA-A*01:59",
                "HLA-A*01:60",
                "HLA-A*01:61",
                "HLA-A*01:62",
                "HLA-A*01:63",
                "HLA-A*01:64",
                "HLA-A*01:65",
                "HLA-A*01:66",
                "HLA-A*02:01",
                "HLA-A*02:02",
                "HLA-A*02:03",
                "HLA-A*02:04",
                "HLA-A*02:05",
                "HLA-A*02:06",
                "HLA-A*02:07",
                "HLA-A*02:08",
                "HLA-A*02:09",
                "HLA-A*02:10",
                "HLA-A*02:101",
                "HLA-A*02:102",
                "HLA-A*02:103",
                "HLA-A*02:104",
                "HLA-A*02:105",
                "HLA-A*02:106",
                "HLA-A*02:107",
                "HLA-A*02:108",
                "HLA-A*02:109",
                "HLA-A*02:11",
                "HLA-A*02:110",
                "HLA-A*02:111",
                "HLA-A*02:112",
                "HLA-A*02:114",
                "HLA-A*02:115",
                "HLA-A*02:116",
                "HLA-A*02:117",
                "HLA-A*02:118",
                "HLA-A*02:119",
                "HLA-A*02:12",
                "HLA-A*02:120",
                "HLA-A*02:121",
                "HLA-A*02:122",
                "HLA-A*02:123",
                "HLA-A*02:124",
                "HLA-A*02:126",
                "HLA-A*02:127",
                "HLA-A*02:128",
                "HLA-A*02:129",
                "HLA-A*02:13",
                "HLA-A*02:130",
                "HLA-A*02:131",
                "HLA-A*02:132",
                "HLA-A*02:133",
                "HLA-A*02:134",
                "HLA-A*02:135",
                "HLA-A*02:136",
                "HLA-A*02:137",
                "HLA-A*02:138",
                "HLA-A*02:139",
                "HLA-A*02:14",
                "HLA-A*02:140",
                "HLA-A*02:141",
                "HLA-A*02:142",
                "HLA-A*02:143",
                "HLA-A*02:144",
                "HLA-A*02:145",
                "HLA-A*02:146",
                "HLA-A*02:147",
                "HLA-A*02:148",
                "HLA-A*02:149",
                "HLA-A*02:150",
                "HLA-A*02:151",
                "HLA-A*02:152",
                "HLA-A*02:153",
                "HLA-A*02:154",
                "HLA-A*02:155",
                "HLA-A*02:156",
                "HLA-A*02:157",
                "HLA-A*02:158",
                "HLA-A*02:159",
                "HLA-A*02:16",
                "HLA-A*02:160",
                "HLA-A*02:161",
                "HLA-A*02:162",
                "HLA-A*02:163",
                "HLA-A*02:164",
                "HLA-A*02:165",
                "HLA-A*02:166",
                "HLA-A*02:167",
                "HLA-A*02:168",
                "HLA-A*02:169",
                "HLA-A*02:17",
                "HLA-A*02:170",
                "HLA-A*02:171",
                "HLA-A*02:172",
                "HLA-A*02:173",
                "HLA-A*02:174",
                "HLA-A*02:175",
                "HLA-A*02:176",
                "HLA-A*02:177",
                "HLA-A*02:178",
                "HLA-A*02:179",
                "HLA-A*02:18",
                "HLA-A*02:180",
                "HLA-A*02:181",
                "HLA-A*02:182",
                "HLA-A*02:183",
                "HLA-A*02:184",
                "HLA-A*02:185",
                "HLA-A*02:186",
                "HLA-A*02:187",
                "HLA-A*02:188",
                "HLA-A*02:189",
                "HLA-A*02:19",
                "HLA-A*02:190",
                "HLA-A*02:191",
                "HLA-A*02:192",
                "HLA-A*02:193",
                "HLA-A*02:194",
                "HLA-A*02:195",
                "HLA-A*02:196",
                "HLA-A*02:197",
                "HLA-A*02:198",
                "HLA-A*02:199",
                "HLA-A*02:20",
                "HLA-A*02:200",
                "HLA-A*02:201",
                "HLA-A*02:202",
                "HLA-A*02:203",
                "HLA-A*02:204",
                "HLA-A*02:205",
                "HLA-A*02:206",
                "HLA-A*02:207",
                "HLA-A*02:208",
                "HLA-A*02:209",
                "HLA-A*02:21",
                "HLA-A*02:210",
                "HLA-A*02:211",
                "HLA-A*02:212",
                "HLA-A*02:213",
                "HLA-A*02:214",
                "HLA-A*02:215",
                "HLA-A*02:216",
                "HLA-A*02:217",
                "HLA-A*02:218",
                "HLA-A*02:219",
                "HLA-A*02:22",
                "HLA-A*02:220",
                "HLA-A*02:221",
                "HLA-A*02:224",
                "HLA-A*02:228",
                "HLA-A*02:229",
                "HLA-A*02:230",
                "HLA-A*02:231",
                "HLA-A*02:232",
                "HLA-A*02:233",
                "HLA-A*02:234",
                "HLA-A*02:235",
                "HLA-A*02:236",
                "HLA-A*02:237",
                "HLA-A*02:238",
                "HLA-A*02:239",
                "HLA-A*02:24",
                "HLA-A*02:240",
                "HLA-A*02:241",
                "HLA-A*02:242",
                "HLA-A*02:243",
                "HLA-A*02:244",
                "HLA-A*02:245",
                "HLA-A*02:246",
                "HLA-A*02:247",
                "HLA-A*02:248",
                "HLA-A*02:249",
                "HLA-A*02:25",
                "HLA-A*02:251",
                "HLA-A*02:252",
                "HLA-A*02:253",
                "HLA-A*02:254",
                "HLA-A*02:255",
                "HLA-A*02:256",
                "HLA-A*02:257",
                "HLA-A*02:258",
                "HLA-A*02:259",
                "HLA-A*02:26",
                "HLA-A*02:260",
                "HLA-A*02:261",
                "HLA-A*02:262",
                "HLA-A*02:263",
                "HLA-A*02:264",
                "HLA-A*02:265",
                "HLA-A*02:266",
                "HLA-A*02:27",
                "HLA-A*02:28",
                "HLA-A*02:29",
                "HLA-A*02:30",
                "HLA-A*02:31",
                "HLA-A*02:33",
                "HLA-A*02:34",
                "HLA-A*02:35",
                "HLA-A*02:36",
                "HLA-A*02:37",
                "HLA-A*02:38",
                "HLA-A*02:39",
                "HLA-A*02:40",
                "HLA-A*02:41",
                "HLA-A*02:42",
                "HLA-A*02:44",
                "HLA-A*02:45",
                "HLA-A*02:46",
                "HLA-A*02:47",
                "HLA-A*02:48",
                "HLA-A*02:49",
                "HLA-A*02:50",
                "HLA-A*02:51",
                "HLA-A*02:52",
                "HLA-A*02:54",
                "HLA-A*02:55",
                "HLA-A*02:56",
                "HLA-A*02:57",
                "HLA-A*02:58",
                "HLA-A*02:59",
                "HLA-A*02:60",
                "HLA-A*02:61",
                "HLA-A*02:62",
                "HLA-A*02:63",
                "HLA-A*02:64",
                "HLA-A*02:65",
                "HLA-A*02:66",
                "HLA-A*02:67",
                "HLA-A*02:68",
                "HLA-A*02:69",
                "HLA-A*02:70",
                "HLA-A*02:71",
                "HLA-A*02:72",
                "HLA-A*02:73",
                "HLA-A*02:74",
                "HLA-A*02:75",
                "HLA-A*02:76",
                "HLA-A*02:77",
                "HLA-A*02:78",
                "HLA-A*02:79",
                "HLA-A*02:80",
                "HLA-A*02:81",
                "HLA-A*02:84",
                "HLA-A*02:85",
                "HLA-A*02:86",
                "HLA-A*02:87",
                "HLA-A*02:89",
                "HLA-A*02:90",
                "HLA-A*02:91",
                "HLA-A*02:92",
                "HLA-A*02:93",
                "HLA-A*02:95",
                "HLA-A*02:96",
                "HLA-A*02:97",
                "HLA-A*02:99",
                "HLA-A*03:01",
                "HLA-A*03:02",
                "HLA-A*03:04",
                "HLA-A*03:05",
                "HLA-A*03:06",
                "HLA-A*03:07",
                "HLA-A*03:08",
                "HLA-A*03:09",
                "HLA-A*03:10",
                "HLA-A*03:12",
                "HLA-A*03:13",
                "HLA-A*03:14",
                "HLA-A*03:15",
                "HLA-A*03:16",
                "HLA-A*03:17",
                "HLA-A*03:18",
                "HLA-A*03:19",
                "HLA-A*03:20",
                "HLA-A*03:22",
                "HLA-A*03:23",
                "HLA-A*03:24",
                "HLA-A*03:25",
                "HLA-A*03:26",
                "HLA-A*03:27",
                "HLA-A*03:28",
                "HLA-A*03:29",
                "HLA-A*03:30",
                "HLA-A*03:31",
                "HLA-A*03:32",
                "HLA-A*03:33",
                "HLA-A*03:34",
                "HLA-A*03:35",
                "HLA-A*03:37",
                "HLA-A*03:38",
                "HLA-A*03:39",
                "HLA-A*03:40",
                "HLA-A*03:41",
                "HLA-A*03:42",
                "HLA-A*03:43",
                "HLA-A*03:44",
                "HLA-A*03:45",
                "HLA-A*03:46",
                "HLA-A*03:47",
                "HLA-A*03:48",
                "HLA-A*03:49",
                "HLA-A*03:50",
                "HLA-A*03:51",
                "HLA-A*03:52",
                "HLA-A*03:53",
                "HLA-A*03:54",
                "HLA-A*03:55",
                "HLA-A*03:56",
                "HLA-A*03:57",
                "HLA-A*03:58",
                "HLA-A*03:59",
                "HLA-A*03:60",
                "HLA-A*03:61",
                "HLA-A*03:62",
                "HLA-A*03:63",
                "HLA-A*03:64",
                "HLA-A*03:65",
                "HLA-A*03:66",
                "HLA-A*03:67",
                "HLA-A*03:70",
                "HLA-A*03:71",
                "HLA-A*03:72",
                "HLA-A*03:73",
                "HLA-A*03:74",
                "HLA-A*03:75",
                "HLA-A*03:76",
                "HLA-A*03:77",
                "HLA-A*03:78",
                "HLA-A*03:79",
                "HLA-A*03:80",
                "HLA-A*03:81",
                "HLA-A*03:82",
                "HLA-A*11:01",
                "HLA-A*11:02",
                "HLA-A*11:03",
                "HLA-A*11:04",
                "HLA-A*11:05",
                "HLA-A*11:06",
                "HLA-A*11:07",
                "HLA-A*11:08",
                "HLA-A*11:09",
                "HLA-A*11:10",
                "HLA-A*11:11",
                "HLA-A*11:12",
                "HLA-A*11:13",
                "HLA-A*11:14",
                "HLA-A*11:15",
                "HLA-A*11:16",
                "HLA-A*11:17",
                "HLA-A*11:18",
                "HLA-A*11:19",
                "HLA-A*11:20",
                "HLA-A*11:22",
                "HLA-A*11:23",
                "HLA-A*11:24",
                "HLA-A*11:25",
                "HLA-A*11:26",
                "HLA-A*11:27",
                "HLA-A*11:29",
                "HLA-A*11:30",
                "HLA-A*11:31",
                "HLA-A*11:32",
                "HLA-A*11:33",
                "HLA-A*11:34",
                "HLA-A*11:35",
                "HLA-A*11:36",
                "HLA-A*11:37",
                "HLA-A*11:38",
                "HLA-A*11:39",
                "HLA-A*11:40",
                "HLA-A*11:41",
                "HLA-A*11:42",
                "HLA-A*11:43",
                "HLA-A*11:44",
                "HLA-A*11:45",
                "HLA-A*11:46",
                "HLA-A*11:47",
                "HLA-A*11:48",
                "HLA-A*11:49",
                "HLA-A*11:51",
                "HLA-A*11:53",
                "HLA-A*11:54",
                "HLA-A*11:55",
                "HLA-A*11:56",
                "HLA-A*11:57",
                "HLA-A*11:58",
                "HLA-A*11:59",
                "HLA-A*11:60",
                "HLA-A*11:61",
                "HLA-A*11:62",
                "HLA-A*11:63",
                "HLA-A*11:64",
                "HLA-A*23:01",
                "HLA-A*23:02",
                "HLA-A*23:03",
                "HLA-A*23:04",
                "HLA-A*23:05",
                "HLA-A*23:06",
                "HLA-A*23:09",
                "HLA-A*23:10",
                "HLA-A*23:12",
                "HLA-A*23:13",
                "HLA-A*23:14",
                "HLA-A*23:15",
                "HLA-A*23:16",
                "HLA-A*23:17",
                "HLA-A*23:18",
                "HLA-A*23:20",
                "HLA-A*23:21",
                "HLA-A*23:22",
                "HLA-A*23:23",
                "HLA-A*23:24",
                "HLA-A*23:25",
                "HLA-A*23:26",
                "HLA-A*24:02",
                "HLA-A*24:03",
                "HLA-A*24:04",
                "HLA-A*24:05",
                "HLA-A*24:06",
                "HLA-A*24:07",
                "HLA-A*24:08",
                "HLA-A*24:10",
                "HLA-A*24:100",
                "HLA-A*24:101",
                "HLA-A*24:102",
                "HLA-A*24:103",
                "HLA-A*24:104",
                "HLA-A*24:105",
                "HLA-A*24:106",
                "HLA-A*24:107",
                "HLA-A*24:108",
                "HLA-A*24:109",
                "HLA-A*24:110",
                "HLA-A*24:111",
                "HLA-A*24:112",
                "HLA-A*24:113",
                "HLA-A*24:114",
                "HLA-A*24:115",
                "HLA-A*24:116",
                "HLA-A*24:117",
                "HLA-A*24:118",
                "HLA-A*24:119",
                "HLA-A*24:120",
                "HLA-A*24:121",
                "HLA-A*24:122",
                "HLA-A*24:123",
                "HLA-A*24:124",
                "HLA-A*24:125",
                "HLA-A*24:126",
                "HLA-A*24:127",
                "HLA-A*24:128",
                "HLA-A*24:129",
                "HLA-A*24:13",
                "HLA-A*24:130",
                "HLA-A*24:131",
                "HLA-A*24:133",
                "HLA-A*24:134",
                "HLA-A*24:135",
                "HLA-A*24:136",
                "HLA-A*24:137",
                "HLA-A*24:138",
                "HLA-A*24:139",
                "HLA-A*24:14",
                "HLA-A*24:140",
                "HLA-A*24:141",
                "HLA-A*24:142",
                "HLA-A*24:143",
                "HLA-A*24:144",
                "HLA-A*24:15",
                "HLA-A*24:17",
                "HLA-A*24:18",
                "HLA-A*24:19",
                "HLA-A*24:20",
                "HLA-A*24:21",
                "HLA-A*24:22",
                "HLA-A*24:23",
                "HLA-A*24:24",
                "HLA-A*24:25",
                "HLA-A*24:26",
                "HLA-A*24:27",
                "HLA-A*24:28",
                "HLA-A*24:29",
                "HLA-A*24:30",
                "HLA-A*24:31",
                "HLA-A*24:32",
                "HLA-A*24:33",
                "HLA-A*24:34",
                "HLA-A*24:35",
                "HLA-A*24:37",
                "HLA-A*24:38",
                "HLA-A*24:39",
                "HLA-A*24:41",
                "HLA-A*24:42",
                "HLA-A*24:43",
                "HLA-A*24:44",
                "HLA-A*24:46",
                "HLA-A*24:47",
                "HLA-A*24:49",
                "HLA-A*24:50",
                "HLA-A*24:51",
                "HLA-A*24:52",
                "HLA-A*24:53",
                "HLA-A*24:54",
                "HLA-A*24:55",
                "HLA-A*24:56",
                "HLA-A*24:57",
                "HLA-A*24:58",
                "HLA-A*24:59",
                "HLA-A*24:61",
                "HLA-A*24:62",
                "HLA-A*24:63",
                "HLA-A*24:64",
                "HLA-A*24:66",
                "HLA-A*24:67",
                "HLA-A*24:68",
                "HLA-A*24:69",
                "HLA-A*24:70",
                "HLA-A*24:71",
                "HLA-A*24:72",
                "HLA-A*24:73",
                "HLA-A*24:74",
                "HLA-A*24:75",
                "HLA-A*24:76",
                "HLA-A*24:77",
                "HLA-A*24:78",
                "HLA-A*24:79",
                "HLA-A*24:80",
                "HLA-A*24:81",
                "HLA-A*24:82",
                "HLA-A*24:85",
                "HLA-A*24:87",
                "HLA-A*24:88",
                "HLA-A*24:89",
                "HLA-A*24:91",
                "HLA-A*24:92",
                "HLA-A*24:93",
                "HLA-A*24:94",
                "HLA-A*24:95",
                "HLA-A*24:96",
                "HLA-A*24:97",
                "HLA-A*24:98",
                "HLA-A*24:99",
                "HLA-A*25:01",
                "HLA-A*25:02",
                "HLA-A*25:03",
                "HLA-A*25:04",
                "HLA-A*25:05",
                "HLA-A*25:06",
                "HLA-A*25:07",
                "HLA-A*25:08",
                "HLA-A*25:09",
                "HLA-A*25:10",
                "HLA-A*25:11",
                "HLA-A*25:13",
                "HLA-A*26:01",
                "HLA-A*26:02",
                "HLA-A*26:03",
                "HLA-A*26:04",
                "HLA-A*26:05",
                "HLA-A*26:06",
                "HLA-A*26:07",
                "HLA-A*26:08",
                "HLA-A*26:09",
                "HLA-A*26:10",
                "HLA-A*26:12",
                "HLA-A*26:13",
                "HLA-A*26:14",
                "HLA-A*26:15",
                "HLA-A*26:16",
                "HLA-A*26:17",
                "HLA-A*26:18",
                "HLA-A*26:19",
                "HLA-A*26:20",
                "HLA-A*26:21",
                "HLA-A*26:22",
                "HLA-A*26:23",
                "HLA-A*26:24",
                "HLA-A*26:26",
                "HLA-A*26:27",
                "HLA-A*26:28",
                "HLA-A*26:29",
                "HLA-A*26:30",
                "HLA-A*26:31",
                "HLA-A*26:32",
                "HLA-A*26:33",
                "HLA-A*26:34",
                "HLA-A*26:35",
                "HLA-A*26:36",
                "HLA-A*26:37",
                "HLA-A*26:38",
                "HLA-A*26:39",
                "HLA-A*26:40",
                "HLA-A*26:41",
                "HLA-A*26:42",
                "HLA-A*26:43",
                "HLA-A*26:45",
                "HLA-A*26:46",
                "HLA-A*26:47",
                "HLA-A*26:48",
                "HLA-A*26:49",
                "HLA-A*26:50",
                "HLA-A*29:01",
                "HLA-A*29:02",
                "HLA-A*29:03",
                "HLA-A*29:04",
                "HLA-A*29:05",
                "HLA-A*29:06",
                "HLA-A*29:07",
                "HLA-A*29:09",
                "HLA-A*29:10",
                "HLA-A*29:11",
                "HLA-A*29:12",
                "HLA-A*29:13",
                "HLA-A*29:14",
                "HLA-A*29:15",
                "HLA-A*29:16",
                "HLA-A*29:17",
                "HLA-A*29:18",
                "HLA-A*29:19",
                "HLA-A*29:20",
                "HLA-A*29:21",
                "HLA-A*29:22",
                "HLA-A*30:01",
                "HLA-A*30:02",
                "HLA-A*30:03",
                "HLA-A*30:04",
                "HLA-A*30:06",
                "HLA-A*30:07",
                "HLA-A*30:08",
                "HLA-A*30:09",
                "HLA-A*30:10",
                "HLA-A*30:11",
                "HLA-A*30:12",
                "HLA-A*30:13",
                "HLA-A*30:15",
                "HLA-A*30:16",
                "HLA-A*30:17",
                "HLA-A*30:18",
                "HLA-A*30:19",
                "HLA-A*30:20",
                "HLA-A*30:22",
                "HLA-A*30:23",
                "HLA-A*30:24",
                "HLA-A*30:25",
                "HLA-A*30:26",
                "HLA-A*30:28",
                "HLA-A*30:29",
                "HLA-A*30:30",
                "HLA-A*30:31",
                "HLA-A*30:32",
                "HLA-A*30:33",
                "HLA-A*30:34",
                "HLA-A*30:35",
                "HLA-A*30:36",
                "HLA-A*30:37",
                "HLA-A*30:38",
                "HLA-A*30:39",
                "HLA-A*30:40",
                "HLA-A*30:41",
                "HLA-A*31:01",
                "HLA-A*31:02",
                "HLA-A*31:03",
                "HLA-A*31:04",
                "HLA-A*31:05",
                "HLA-A*31:06",
                "HLA-A*31:07",
                "HLA-A*31:08",
                "HLA-A*31:09",
                "HLA-A*31:10",
                "HLA-A*31:11",
                "HLA-A*31:12",
                "HLA-A*31:13",
                "HLA-A*31:15",
                "HLA-A*31:16",
                "HLA-A*31:17",
                "HLA-A*31:18",
                "HLA-A*31:19",
                "HLA-A*31:20",
                "HLA-A*31:21",
                "HLA-A*31:22",
                "HLA-A*31:23",
                "HLA-A*31:24",
                "HLA-A*31:25",
                "HLA-A*31:26",
                "HLA-A*31:27",
                "HLA-A*31:28",
                "HLA-A*31:29",
                "HLA-A*31:30",
                "HLA-A*31:31",
                "HLA-A*31:32",
                "HLA-A*31:33",
                "HLA-A*31:34",
                "HLA-A*31:35",
                "HLA-A*31:36",
                "HLA-A*31:37",
                "HLA-A*32:01",
                "HLA-A*32:02",
                "HLA-A*32:03",
                "HLA-A*32:04",
                "HLA-A*32:05",
                "HLA-A*32:06",
                "HLA-A*32:07",
                "HLA-A*32:08",
                "HLA-A*32:09",
                "HLA-A*32:10",
                "HLA-A*32:12",
                "HLA-A*32:13",
                "HLA-A*32:14",
                "HLA-A*32:15",
                "HLA-A*32:16",
                "HLA-A*32:17",
                "HLA-A*32:18",
                "HLA-A*32:20",
                "HLA-A*32:21",
                "HLA-A*32:22",
                "HLA-A*32:23",
                "HLA-A*32:24",
                "HLA-A*32:25",
                "HLA-A*33:01",
                "HLA-A*33:03",
                "HLA-A*33:04",
                "HLA-A*33:05",
                "HLA-A*33:06",
                "HLA-A*33:07",
                "HLA-A*33:08",
                "HLA-A*33:09",
                "HLA-A*33:10",
                "HLA-A*33:11",
                "HLA-A*33:12",
                "HLA-A*33:13",
                "HLA-A*33:14",
                "HLA-A*33:15",
                "HLA-A*33:16",
                "HLA-A*33:17",
                "HLA-A*33:18",
                "HLA-A*33:19",
                "HLA-A*33:20",
                "HLA-A*33:21",
                "HLA-A*33:22",
                "HLA-A*33:23",
                "HLA-A*33:24",
                "HLA-A*33:25",
                "HLA-A*33:26",
                "HLA-A*33:27",
                "HLA-A*33:28",
                "HLA-A*33:29",
                "HLA-A*33:30",
                "HLA-A*33:31",
                "HLA-A*34:01",
                "HLA-A*34:02",
                "HLA-A*34:03",
                "HLA-A*34:04",
                "HLA-A*34:05",
                "HLA-A*34:06",
                "HLA-A*34:07",
                "HLA-A*34:08",
                "HLA-A*36:01",
                "HLA-A*36:02",
                "HLA-A*36:03",
                "HLA-A*36:04",
                "HLA-A*36:05",
                "HLA-A*43:01",
                "HLA-A*66:01",
                "HLA-A*66:02",
                "HLA-A*66:03",
                "HLA-A*66:04",
                "HLA-A*66:05",
                "HLA-A*66:06",
                "HLA-A*66:07",
                "HLA-A*66:08",
                "HLA-A*66:09",
                "HLA-A*66:10",
                "HLA-A*66:11",
                "HLA-A*66:12",
                "HLA-A*66:13",
                "HLA-A*66:14",
                "HLA-A*66:15",
                "HLA-A*68:01",
                "HLA-A*68:02",
                "HLA-A*68:03",
                "HLA-A*68:04",
                "HLA-A*68:05",
                "HLA-A*68:06",
                "HLA-A*68:07",
                "HLA-A*68:08",
                "HLA-A*68:09",
                "HLA-A*68:10",
                "HLA-A*68:12",
                "HLA-A*68:13",
                "HLA-A*68:14",
                "HLA-A*68:15",
                "HLA-A*68:16",
                "HLA-A*68:17",
                "HLA-A*68:19",
                "HLA-A*68:20",
                "HLA-A*68:21",
                "HLA-A*68:22",
                "HLA-A*68:23",
                "HLA-A*68:24",
                "HLA-A*68:25",
                "HLA-A*68:26",
                "HLA-A*68:27",
                "HLA-A*68:28",
                "HLA-A*68:29",
                "HLA-A*68:30",
                "HLA-A*68:31",
                "HLA-A*68:32",
                "HLA-A*68:33",
                "HLA-A*68:34",
                "HLA-A*68:35",
                "HLA-A*68:36",
                "HLA-A*68:37",
                "HLA-A*68:38",
                "HLA-A*68:39",
                "HLA-A*68:40",
                "HLA-A*68:41",
                "HLA-A*68:42",
                "HLA-A*68:43",
                "HLA-A*68:44",
                "HLA-A*68:45",
                "HLA-A*68:46",
                "HLA-A*68:47",
                "HLA-A*68:48",
                "HLA-A*68:50",
                "HLA-A*68:51",
                "HLA-A*68:52",
                "HLA-A*68:53",
                "HLA-A*68:54",
                "HLA-A*69:01",
                "HLA-A*74:01",
                "HLA-A*74:02",
                "HLA-A*74:03",
                "HLA-A*74:04",
                "HLA-A*74:05",
                "HLA-A*74:06",
                "HLA-A*74:07",
                "HLA-A*74:08",
                "HLA-A*74:09",
                "HLA-A*74:10",
                "HLA-A*74:11",
                "HLA-A*74:13",
                "HLA-A*80:01",
                "HLA-A*80:02",
                "HLA-B*07:02",
                "HLA-B*07:03",
                "HLA-B*07:04",
                "HLA-B*07:05",
                "HLA-B*07:06",
                "HLA-B*07:07",
                "HLA-B*07:08",
                "HLA-B*07:09",
                "HLA-B*07:10",
                "HLA-B*07:100",
                "HLA-B*07:101",
                "HLA-B*07:102",
                "HLA-B*07:103",
                "HLA-B*07:104",
                "HLA-B*07:105",
                "HLA-B*07:106",
                "HLA-B*07:107",
                "HLA-B*07:108",
                "HLA-B*07:109",
                "HLA-B*07:11",
                "HLA-B*07:110",
                "HLA-B*07:112",
                "HLA-B*07:113",
                "HLA-B*07:114",
                "HLA-B*07:115",
                "HLA-B*07:12",
                "HLA-B*07:13",
                "HLA-B*07:14",
                "HLA-B*07:15",
                "HLA-B*07:16",
                "HLA-B*07:17",
                "HLA-B*07:18",
                "HLA-B*07:19",
                "HLA-B*07:20",
                "HLA-B*07:21",
                "HLA-B*07:22",
                "HLA-B*07:23",
                "HLA-B*07:24",
                "HLA-B*07:25",
                "HLA-B*07:26",
                "HLA-B*07:27",
                "HLA-B*07:28",
                "HLA-B*07:29",
                "HLA-B*07:30",
                "HLA-B*07:31",
                "HLA-B*07:32",
                "HLA-B*07:33",
                "HLA-B*07:34",
                "HLA-B*07:35",
                "HLA-B*07:36",
                "HLA-B*07:37",
                "HLA-B*07:38",
                "HLA-B*07:39",
                "HLA-B*07:40",
                "HLA-B*07:41",
                "HLA-B*07:42",
                "HLA-B*07:43",
                "HLA-B*07:44",
                "HLA-B*07:45",
                "HLA-B*07:46",
                "HLA-B*07:47",
                "HLA-B*07:48",
                "HLA-B*07:50",
                "HLA-B*07:51",
                "HLA-B*07:52",
                "HLA-B*07:53",
                "HLA-B*07:54",
                "HLA-B*07:55",
                "HLA-B*07:56",
                "HLA-B*07:57",
                "HLA-B*07:58",
                "HLA-B*07:59",
                "HLA-B*07:60",
                "HLA-B*07:61",
                "HLA-B*07:62",
                "HLA-B*07:63",
                "HLA-B*07:64",
                "HLA-B*07:65",
                "HLA-B*07:66",
                "HLA-B*07:68",
                "HLA-B*07:69",
                "HLA-B*07:70",
                "HLA-B*07:71",
                "HLA-B*07:72",
                "HLA-B*07:73",
                "HLA-B*07:74",
                "HLA-B*07:75",
                "HLA-B*07:76",
                "HLA-B*07:77",
                "HLA-B*07:78",
                "HLA-B*07:79",
                "HLA-B*07:80",
                "HLA-B*07:81",
                "HLA-B*07:82",
                "HLA-B*07:83",
                "HLA-B*07:84",
                "HLA-B*07:85",
                "HLA-B*07:86",
                "HLA-B*07:87",
                "HLA-B*07:88",
                "HLA-B*07:89",
                "HLA-B*07:90",
                "HLA-B*07:91",
                "HLA-B*07:92",
                "HLA-B*07:93",
                "HLA-B*07:94",
                "HLA-B*07:95",
                "HLA-B*07:96",
                "HLA-B*07:97",
                "HLA-B*07:98",
                "HLA-B*07:99",
                "HLA-B*08:01",
                "HLA-B*08:02",
                "HLA-B*08:03",
                "HLA-B*08:04",
                "HLA-B*08:05",
                "HLA-B*08:07",
                "HLA-B*08:09",
                "HLA-B*08:10",
                "HLA-B*08:11",
                "HLA-B*08:12",
                "HLA-B*08:13",
                "HLA-B*08:14",
                "HLA-B*08:15",
                "HLA-B*08:16",
                "HLA-B*08:17",
                "HLA-B*08:18",
                "HLA-B*08:20",
                "HLA-B*08:21",
                "HLA-B*08:22",
                "HLA-B*08:23",
                "HLA-B*08:24",
                "HLA-B*08:25",
                "HLA-B*08:26",
                "HLA-B*08:27",
                "HLA-B*08:28",
                "HLA-B*08:29",
                "HLA-B*08:31",
                "HLA-B*08:32",
                "HLA-B*08:33",
                "HLA-B*08:34",
                "HLA-B*08:35",
                "HLA-B*08:36",
                "HLA-B*08:37",
                "HLA-B*08:38",
                "HLA-B*08:39",
                "HLA-B*08:40",
                "HLA-B*08:41",
                "HLA-B*08:42",
                "HLA-B*08:43",
                "HLA-B*08:44",
                "HLA-B*08:45",
                "HLA-B*08:46",
                "HLA-B*08:47",
                "HLA-B*08:48",
                "HLA-B*08:49",
                "HLA-B*08:50",
                "HLA-B*08:51",
                "HLA-B*08:52",
                "HLA-B*08:53",
                "HLA-B*08:54",
                "HLA-B*08:55",
                "HLA-B*08:56",
                "HLA-B*08:57",
                "HLA-B*08:58",
                "HLA-B*08:59",
                "HLA-B*08:60",
                "HLA-B*08:61",
                "HLA-B*08:62",
                "HLA-B*13:01",
                "HLA-B*13:02",
                "HLA-B*13:03",
                "HLA-B*13:04",
                "HLA-B*13:06",
                "HLA-B*13:09",
                "HLA-B*13:10",
                "HLA-B*13:11",
                "HLA-B*13:12",
                "HLA-B*13:13",
                "HLA-B*13:14",
                "HLA-B*13:15",
                "HLA-B*13:16",
                "HLA-B*13:17",
                "HLA-B*13:18",
                "HLA-B*13:19",
                "HLA-B*13:20",
                "HLA-B*13:21",
                "HLA-B*13:22",
                "HLA-B*13:23",
                "HLA-B*13:25",
                "HLA-B*13:26",
                "HLA-B*13:27",
                "HLA-B*13:28",
                "HLA-B*13:29",
                "HLA-B*13:30",
                "HLA-B*13:31",
                "HLA-B*13:32",
                "HLA-B*13:33",
                "HLA-B*13:34",
                "HLA-B*13:35",
                "HLA-B*13:36",
                "HLA-B*13:37",
                "HLA-B*13:38",
                "HLA-B*13:39",
                "HLA-B*14:01",
                "HLA-B*14:02",
                "HLA-B*14:03",
                "HLA-B*14:04",
                "HLA-B*14:05",
                "HLA-B*14:06",
                "HLA-B*14:08",
                "HLA-B*14:09",
                "HLA-B*14:10",
                "HLA-B*14:11",
                "HLA-B*14:12",
                "HLA-B*14:13",
                "HLA-B*14:14",
                "HLA-B*14:15",
                "HLA-B*14:16",
                "HLA-B*14:17",
                "HLA-B*14:18",
                "HLA-B*15:01",
                "HLA-B*15:02",
                "HLA-B*15:03",
                "HLA-B*15:04",
                "HLA-B*15:05",
                "HLA-B*15:06",
                "HLA-B*15:07",
                "HLA-B*15:08",
                "HLA-B*15:09",
                "HLA-B*15:10",
                "HLA-B*15:101",
                "HLA-B*15:102",
                "HLA-B*15:103",
                "HLA-B*15:104",
                "HLA-B*15:105",
                "HLA-B*15:106",
                "HLA-B*15:107",
                "HLA-B*15:108",
                "HLA-B*15:109",
                "HLA-B*15:11",
                "HLA-B*15:110",
                "HLA-B*15:112",
                "HLA-B*15:113",
                "HLA-B*15:114",
                "HLA-B*15:115",
                "HLA-B*15:116",
                "HLA-B*15:117",
                "HLA-B*15:118",
                "HLA-B*15:119",
                "HLA-B*15:12",
                "HLA-B*15:120",
                "HLA-B*15:121",
                "HLA-B*15:122",
                "HLA-B*15:123",
                "HLA-B*15:124",
                "HLA-B*15:125",
                "HLA-B*15:126",
                "HLA-B*15:127",
                "HLA-B*15:128",
                "HLA-B*15:129",
                "HLA-B*15:13",
                "HLA-B*15:131",
                "HLA-B*15:132",
                "HLA-B*15:133",
                "HLA-B*15:134",
                "HLA-B*15:135",
                "HLA-B*15:136",
                "HLA-B*15:137",
                "HLA-B*15:138",
                "HLA-B*15:139",
                "HLA-B*15:14",
                "HLA-B*15:140",
                "HLA-B*15:141",
                "HLA-B*15:142",
                "HLA-B*15:143",
                "HLA-B*15:144",
                "HLA-B*15:145",
                "HLA-B*15:146",
                "HLA-B*15:147",
                "HLA-B*15:148",
                "HLA-B*15:15",
                "HLA-B*15:150",
                "HLA-B*15:151",
                "HLA-B*15:152",
                "HLA-B*15:153",
                "HLA-B*15:154",
                "HLA-B*15:155",
                "HLA-B*15:156",
                "HLA-B*15:157",
                "HLA-B*15:158",
                "HLA-B*15:159",
                "HLA-B*15:16",
                "HLA-B*15:160",
                "HLA-B*15:161",
                "HLA-B*15:162",
                "HLA-B*15:163",
                "HLA-B*15:164",
                "HLA-B*15:165",
                "HLA-B*15:166",
                "HLA-B*15:167",
                "HLA-B*15:168",
                "HLA-B*15:169",
                "HLA-B*15:17",
                "HLA-B*15:170",
                "HLA-B*15:171",
                "HLA-B*15:172",
                "HLA-B*15:173",
                "HLA-B*15:174",
                "HLA-B*15:175",
                "HLA-B*15:176",
                "HLA-B*15:177",
                "HLA-B*15:178",
                "HLA-B*15:179",
                "HLA-B*15:18",
                "HLA-B*15:180",
                "HLA-B*15:183",
                "HLA-B*15:184",
                "HLA-B*15:185",
                "HLA-B*15:186",
                "HLA-B*15:187",
                "HLA-B*15:188",
                "HLA-B*15:189",
                "HLA-B*15:19",
                "HLA-B*15:191",
                "HLA-B*15:192",
                "HLA-B*15:193",
                "HLA-B*15:194",
                "HLA-B*15:195",
                "HLA-B*15:196",
                "HLA-B*15:197",
                "HLA-B*15:198",
                "HLA-B*15:199",
                "HLA-B*15:20",
                "HLA-B*15:200",
                "HLA-B*15:201",
                "HLA-B*15:202",
                "HLA-B*15:21",
                "HLA-B*15:23",
                "HLA-B*15:24",
                "HLA-B*15:25",
                "HLA-B*15:27",
                "HLA-B*15:28",
                "HLA-B*15:29",
                "HLA-B*15:30",
                "HLA-B*15:31",
                "HLA-B*15:32",
                "HLA-B*15:33",
                "HLA-B*15:34",
                "HLA-B*15:35",
                "HLA-B*15:36",
                "HLA-B*15:37",
                "HLA-B*15:38",
                "HLA-B*15:39",
                "HLA-B*15:40",
                "HLA-B*15:42",
                "HLA-B*15:43",
                "HLA-B*15:44",
                "HLA-B*15:45",
                "HLA-B*15:46",
                "HLA-B*15:47",
                "HLA-B*15:48",
                "HLA-B*15:49",
                "HLA-B*15:50",
                "HLA-B*15:51",
                "HLA-B*15:52",
                "HLA-B*15:53",
                "HLA-B*15:54",
                "HLA-B*15:55",
                "HLA-B*15:56",
                "HLA-B*15:57",
                "HLA-B*15:58",
                "HLA-B*15:60",
                "HLA-B*15:61",
                "HLA-B*15:62",
                "HLA-B*15:63",
                "HLA-B*15:64",
                "HLA-B*15:65",
                "HLA-B*15:66",
                "HLA-B*15:67",
                "HLA-B*15:68",
                "HLA-B*15:69",
                "HLA-B*15:70",
                "HLA-B*15:71",
                "HLA-B*15:72",
                "HLA-B*15:73",
                "HLA-B*15:74",
                "HLA-B*15:75",
                "HLA-B*15:76",
                "HLA-B*15:77",
                "HLA-B*15:78",
                "HLA-B*15:80",
                "HLA-B*15:81",
                "HLA-B*15:82",
                "HLA-B*15:83",
                "HLA-B*15:84",
                "HLA-B*15:85",
                "HLA-B*15:86",
                "HLA-B*15:87",
                "HLA-B*15:88",
                "HLA-B*15:89",
                "HLA-B*15:90",
                "HLA-B*15:91",
                "HLA-B*15:92",
                "HLA-B*15:93",
                "HLA-B*15:95",
                "HLA-B*15:96",
                "HLA-B*15:97",
                "HLA-B*15:98",
                "HLA-B*15:99",
                "HLA-B*18:01",
                "HLA-B*18:02",
                "HLA-B*18:03",
                "HLA-B*18:04",
                "HLA-B*18:05",
                "HLA-B*18:06",
                "HLA-B*18:07",
                "HLA-B*18:08",
                "HLA-B*18:09",
                "HLA-B*18:10",
                "HLA-B*18:11",
                "HLA-B*18:12",
                "HLA-B*18:13",
                "HLA-B*18:14",
                "HLA-B*18:15",
                "HLA-B*18:18",
                "HLA-B*18:19",
                "HLA-B*18:20",
                "HLA-B*18:21",
                "HLA-B*18:22",
                "HLA-B*18:24",
                "HLA-B*18:25",
                "HLA-B*18:26",
                "HLA-B*18:27",
                "HLA-B*18:28",
                "HLA-B*18:29",
                "HLA-B*18:30",
                "HLA-B*18:31",
                "HLA-B*18:32",
                "HLA-B*18:33",
                "HLA-B*18:34",
                "HLA-B*18:35",
                "HLA-B*18:36",
                "HLA-B*18:37",
                "HLA-B*18:38",
                "HLA-B*18:39",
                "HLA-B*18:40",
                "HLA-B*18:41",
                "HLA-B*18:42",
                "HLA-B*18:43",
                "HLA-B*18:44",
                "HLA-B*18:45",
                "HLA-B*18:46",
                "HLA-B*18:47",
                "HLA-B*18:48",
                "HLA-B*18:49",
                "HLA-B*18:50",
                "HLA-B*27:01",
                "HLA-B*27:02",
                "HLA-B*27:03",
                "HLA-B*27:04",
                "HLA-B*27:05",
                "HLA-B*27:06",
                "HLA-B*27:07",
                "HLA-B*27:08",
                "HLA-B*27:09",
                "HLA-B*27:10",
                "HLA-B*27:11",
                "HLA-B*27:12",
                "HLA-B*27:13",
                "HLA-B*27:14",
                "HLA-B*27:15",
                "HLA-B*27:16",
                "HLA-B*27:17",
                "HLA-B*27:18",
                "HLA-B*27:19",
                "HLA-B*27:20",
                "HLA-B*27:21",
                "HLA-B*27:23",
                "HLA-B*27:24",
                "HLA-B*27:25",
                "HLA-B*27:26",
                "HLA-B*27:27",
                "HLA-B*27:28",
                "HLA-B*27:29",
                "HLA-B*27:30",
                "HLA-B*27:31",
                "HLA-B*27:32",
                "HLA-B*27:33",
                "HLA-B*27:34",
                "HLA-B*27:35",
                "HLA-B*27:36",
                "HLA-B*27:37",
                "HLA-B*27:38",
                "HLA-B*27:39",
                "HLA-B*27:40",
                "HLA-B*27:41",
                "HLA-B*27:42",
                "HLA-B*27:43",
                "HLA-B*27:44",
                "HLA-B*27:45",
                "HLA-B*27:46",
                "HLA-B*27:47",
                "HLA-B*27:48",
                "HLA-B*27:49",
                "HLA-B*27:50",
                "HLA-B*27:51",
                "HLA-B*27:52",
                "HLA-B*27:53",
                "HLA-B*27:54",
                "HLA-B*27:55",
                "HLA-B*27:56",
                "HLA-B*27:57",
                "HLA-B*27:58",
                "HLA-B*27:60",
                "HLA-B*27:61",
                "HLA-B*27:62",
                "HLA-B*27:63",
                "HLA-B*27:67",
                "HLA-B*27:68",
                "HLA-B*27:69",
                "HLA-B*35:01",
                "HLA-B*35:02",
                "HLA-B*35:03",
                "HLA-B*35:04",
                "HLA-B*35:05",
                "HLA-B*35:06",
                "HLA-B*35:07",
                "HLA-B*35:08",
                "HLA-B*35:09",
                "HLA-B*35:10",
                "HLA-B*35:100",
                "HLA-B*35:101",
                "HLA-B*35:102",
                "HLA-B*35:103",
                "HLA-B*35:104",
                "HLA-B*35:105",
                "HLA-B*35:106",
                "HLA-B*35:107",
                "HLA-B*35:108",
                "HLA-B*35:109",
                "HLA-B*35:11",
                "HLA-B*35:110",
                "HLA-B*35:111",
                "HLA-B*35:112",
                "HLA-B*35:113",
                "HLA-B*35:114",
                "HLA-B*35:115",
                "HLA-B*35:116",
                "HLA-B*35:117",
                "HLA-B*35:118",
                "HLA-B*35:119",
                "HLA-B*35:12",
                "HLA-B*35:120",
                "HLA-B*35:121",
                "HLA-B*35:122",
                "HLA-B*35:123",
                "HLA-B*35:124",
                "HLA-B*35:125",
                "HLA-B*35:126",
                "HLA-B*35:127",
                "HLA-B*35:128",
                "HLA-B*35:13",
                "HLA-B*35:131",
                "HLA-B*35:132",
                "HLA-B*35:133",
                "HLA-B*35:135",
                "HLA-B*35:136",
                "HLA-B*35:137",
                "HLA-B*35:138",
                "HLA-B*35:139",
                "HLA-B*35:14",
                "HLA-B*35:140",
                "HLA-B*35:141",
                "HLA-B*35:142",
                "HLA-B*35:143",
                "HLA-B*35:144",
                "HLA-B*35:15",
                "HLA-B*35:16",
                "HLA-B*35:17",
                "HLA-B*35:18",
                "HLA-B*35:19",
                "HLA-B*35:20",
                "HLA-B*35:21",
                "HLA-B*35:22",
                "HLA-B*35:23",
                "HLA-B*35:24",
                "HLA-B*35:25",
                "HLA-B*35:26",
                "HLA-B*35:27",
                "HLA-B*35:28",
                "HLA-B*35:29",
                "HLA-B*35:30",
                "HLA-B*35:31",
                "HLA-B*35:32",
                "HLA-B*35:33",
                "HLA-B*35:34",
                "HLA-B*35:35",
                "HLA-B*35:36",
                "HLA-B*35:37",
                "HLA-B*35:38",
                "HLA-B*35:39",
                "HLA-B*35:41",
                "HLA-B*35:42",
                "HLA-B*35:43",
                "HLA-B*35:44",
                "HLA-B*35:45",
                "HLA-B*35:46",
                "HLA-B*35:47",
                "HLA-B*35:48",
                "HLA-B*35:49",
                "HLA-B*35:50",
                "HLA-B*35:51",
                "HLA-B*35:52",
                "HLA-B*35:54",
                "HLA-B*35:55",
                "HLA-B*35:56",
                "HLA-B*35:57",
                "HLA-B*35:58",
                "HLA-B*35:59",
                "HLA-B*35:60",
                "HLA-B*35:61",
                "HLA-B*35:62",
                "HLA-B*35:63",
                "HLA-B*35:64",
                "HLA-B*35:66",
                "HLA-B*35:67",
                "HLA-B*35:68",
                "HLA-B*35:69",
                "HLA-B*35:70",
                "HLA-B*35:71",
                "HLA-B*35:72",
                "HLA-B*35:74",
                "HLA-B*35:75",
                "HLA-B*35:76",
                "HLA-B*35:77",
                "HLA-B*35:78",
                "HLA-B*35:79",
                "HLA-B*35:80",
                "HLA-B*35:81",
                "HLA-B*35:82",
                "HLA-B*35:83",
                "HLA-B*35:84",
                "HLA-B*35:85",
                "HLA-B*35:86",
                "HLA-B*35:87",
                "HLA-B*35:88",
                "HLA-B*35:89",
                "HLA-B*35:90",
                "HLA-B*35:91",
                "HLA-B*35:92",
                "HLA-B*35:93",
                "HLA-B*35:94",
                "HLA-B*35:95",
                "HLA-B*35:96",
                "HLA-B*35:97",
                "HLA-B*35:98",
                "HLA-B*35:99",
                "HLA-B*37:01",
                "HLA-B*37:02",
                "HLA-B*37:04",
                "HLA-B*37:05",
                "HLA-B*37:06",
                "HLA-B*37:07",
                "HLA-B*37:08",
                "HLA-B*37:09",
                "HLA-B*37:10",
                "HLA-B*37:11",
                "HLA-B*37:12",
                "HLA-B*37:13",
                "HLA-B*37:14",
                "HLA-B*37:15",
                "HLA-B*37:17",
                "HLA-B*37:18",
                "HLA-B*37:19",
                "HLA-B*37:20",
                "HLA-B*37:21",
                "HLA-B*37:22",
                "HLA-B*37:23",
                "HLA-B*38:01",
                "HLA-B*38:02",
                "HLA-B*38:03",
                "HLA-B*38:04",
                "HLA-B*38:05",
                "HLA-B*38:06",
                "HLA-B*38:07",
                "HLA-B*38:08",
                "HLA-B*38:09",
                "HLA-B*38:10",
                "HLA-B*38:11",
                "HLA-B*38:12",
                "HLA-B*38:13",
                "HLA-B*38:14",
                "HLA-B*38:15",
                "HLA-B*38:16",
                "HLA-B*38:17",
                "HLA-B*38:18",
                "HLA-B*38:19",
                "HLA-B*38:20",
                "HLA-B*38:21",
                "HLA-B*38:22",
                "HLA-B*38:23",
                "HLA-B*39:01",
                "HLA-B*39:02",
                "HLA-B*39:03",
                "HLA-B*39:04",
                "HLA-B*39:05",
                "HLA-B*39:06",
                "HLA-B*39:07",
                "HLA-B*39:08",
                "HLA-B*39:09",
                "HLA-B*39:10",
                "HLA-B*39:11",
                "HLA-B*39:12",
                "HLA-B*39:13",
                "HLA-B*39:14",
                "HLA-B*39:15",
                "HLA-B*39:16",
                "HLA-B*39:17",
                "HLA-B*39:18",
                "HLA-B*39:19",
                "HLA-B*39:20",
                "HLA-B*39:22",
                "HLA-B*39:23",
                "HLA-B*39:24",
                "HLA-B*39:26",
                "HLA-B*39:27",
                "HLA-B*39:28",
                "HLA-B*39:29",
                "HLA-B*39:30",
                "HLA-B*39:31",
                "HLA-B*39:32",
                "HLA-B*39:33",
                "HLA-B*39:34",
                "HLA-B*39:35",
                "HLA-B*39:36",
                "HLA-B*39:37",
                "HLA-B*39:39",
                "HLA-B*39:41",
                "HLA-B*39:42",
                "HLA-B*39:43",
                "HLA-B*39:44",
                "HLA-B*39:45",
                "HLA-B*39:46",
                "HLA-B*39:47",
                "HLA-B*39:48",
                "HLA-B*39:49",
                "HLA-B*39:50",
                "HLA-B*39:51",
                "HLA-B*39:52",
                "HLA-B*39:53",
                "HLA-B*39:54",
                "HLA-B*39:55",
                "HLA-B*39:56",
                "HLA-B*39:57",
                "HLA-B*39:58",
                "HLA-B*39:59",
                "HLA-B*39:60",
                "HLA-B*40:01",
                "HLA-B*40:02",
                "HLA-B*40:03",
                "HLA-B*40:04",
                "HLA-B*40:05",
                "HLA-B*40:06",
                "HLA-B*40:07",
                "HLA-B*40:08",
                "HLA-B*40:09",
                "HLA-B*40:10",
                "HLA-B*40:100",
                "HLA-B*40:101",
                "HLA-B*40:102",
                "HLA-B*40:103",
                "HLA-B*40:104",
                "HLA-B*40:105",
                "HLA-B*40:106",
                "HLA-B*40:107",
                "HLA-B*40:108",
                "HLA-B*40:109",
                "HLA-B*40:11",
                "HLA-B*40:110",
                "HLA-B*40:111",
                "HLA-B*40:112",
                "HLA-B*40:113",
                "HLA-B*40:114",
                "HLA-B*40:115",
                "HLA-B*40:116",
                "HLA-B*40:117",
                "HLA-B*40:119",
                "HLA-B*40:12",
                "HLA-B*40:120",
                "HLA-B*40:121",
                "HLA-B*40:122",
                "HLA-B*40:123",
                "HLA-B*40:124",
                "HLA-B*40:125",
                "HLA-B*40:126",
                "HLA-B*40:127",
                "HLA-B*40:128",
                "HLA-B*40:129",
                "HLA-B*40:13",
                "HLA-B*40:130",
                "HLA-B*40:131",
                "HLA-B*40:132",
                "HLA-B*40:134",
                "HLA-B*40:135",
                "HLA-B*40:136",
                "HLA-B*40:137",
                "HLA-B*40:138",
                "HLA-B*40:139",
                "HLA-B*40:14",
                "HLA-B*40:140",
                "HLA-B*40:141",
                "HLA-B*40:143",
                "HLA-B*40:145",
                "HLA-B*40:146",
                "HLA-B*40:147",
                "HLA-B*40:15",
                "HLA-B*40:16",
                "HLA-B*40:18",
                "HLA-B*40:19",
                "HLA-B*40:20",
                "HLA-B*40:21",
                "HLA-B*40:23",
                "HLA-B*40:24",
                "HLA-B*40:25",
                "HLA-B*40:26",
                "HLA-B*40:27",
                "HLA-B*40:28",
                "HLA-B*40:29",
                "HLA-B*40:30",
                "HLA-B*40:31",
                "HLA-B*40:32",
                "HLA-B*40:33",
                "HLA-B*40:34",
                "HLA-B*40:35",
                "HLA-B*40:36",
                "HLA-B*40:37",
                "HLA-B*40:38",
                "HLA-B*40:39",
                "HLA-B*40:40",
                "HLA-B*40:42",
                "HLA-B*40:43",
                "HLA-B*40:44",
                "HLA-B*40:45",
                "HLA-B*40:46",
                "HLA-B*40:47",
                "HLA-B*40:48",
                "HLA-B*40:49",
                "HLA-B*40:50",
                "HLA-B*40:51",
                "HLA-B*40:52",
                "HLA-B*40:53",
                "HLA-B*40:54",
                "HLA-B*40:55",
                "HLA-B*40:56",
                "HLA-B*40:57",
                "HLA-B*40:58",
                "HLA-B*40:59",
                "HLA-B*40:60",
                "HLA-B*40:61",
                "HLA-B*40:62",
                "HLA-B*40:63",
                "HLA-B*40:64",
                "HLA-B*40:65",
                "HLA-B*40:66",
                "HLA-B*40:67",
                "HLA-B*40:68",
                "HLA-B*40:69",
                "HLA-B*40:70",
                "HLA-B*40:71",
                "HLA-B*40:72",
                "HLA-B*40:73",
                "HLA-B*40:74",
                "HLA-B*40:75",
                "HLA-B*40:76",
                "HLA-B*40:77",
                "HLA-B*40:78",
                "HLA-B*40:79",
                "HLA-B*40:80",
                "HLA-B*40:81",
                "HLA-B*40:82",
                "HLA-B*40:83",
                "HLA-B*40:84",
                "HLA-B*40:85",
                "HLA-B*40:86",
                "HLA-B*40:87",
                "HLA-B*40:88",
                "HLA-B*40:89",
                "HLA-B*40:90",
                "HLA-B*40:91",
                "HLA-B*40:92",
                "HLA-B*40:93",
                "HLA-B*40:94",
                "HLA-B*40:95",
                "HLA-B*40:96",
                "HLA-B*40:97",
                "HLA-B*40:98",
                "HLA-B*40:99",
                "HLA-B*41:01",
                "HLA-B*41:02",
                "HLA-B*41:03",
                "HLA-B*41:04",
                "HLA-B*41:05",
                "HLA-B*41:06",
                "HLA-B*41:07",
                "HLA-B*41:08",
                "HLA-B*41:09",
                "HLA-B*41:10",
                "HLA-B*41:11",
                "HLA-B*41:12",
                "HLA-B*42:01",
                "HLA-B*42:02",
                "HLA-B*42:04",
                "HLA-B*42:05",
                "HLA-B*42:06",
                "HLA-B*42:07",
                "HLA-B*42:08",
                "HLA-B*42:09",
                "HLA-B*42:10",
                "HLA-B*42:11",
                "HLA-B*42:12",
                "HLA-B*42:13",
                "HLA-B*42:14",
                "HLA-B*44:02",
                "HLA-B*44:03",
                "HLA-B*44:04",
                "HLA-B*44:05",
                "HLA-B*44:06",
                "HLA-B*44:07",
                "HLA-B*44:08",
                "HLA-B*44:09",
                "HLA-B*44:10",
                "HLA-B*44:100",
                "HLA-B*44:101",
                "HLA-B*44:102",
                "HLA-B*44:103",
                "HLA-B*44:104",
                "HLA-B*44:105",
                "HLA-B*44:106",
                "HLA-B*44:107",
                "HLA-B*44:109",
                "HLA-B*44:11",
                "HLA-B*44:110",
                "HLA-B*44:12",
                "HLA-B*44:13",
                "HLA-B*44:14",
                "HLA-B*44:15",
                "HLA-B*44:16",
                "HLA-B*44:17",
                "HLA-B*44:18",
                "HLA-B*44:20",
                "HLA-B*44:21",
                "HLA-B*44:22",
                "HLA-B*44:24",
                "HLA-B*44:25",
                "HLA-B*44:26",
                "HLA-B*44:27",
                "HLA-B*44:28",
                "HLA-B*44:29",
                "HLA-B*44:30",
                "HLA-B*44:31",
                "HLA-B*44:32",
                "HLA-B*44:33",
                "HLA-B*44:34",
                "HLA-B*44:35",
                "HLA-B*44:36",
                "HLA-B*44:37",
                "HLA-B*44:38",
                "HLA-B*44:39",
                "HLA-B*44:40",
                "HLA-B*44:41",
                "HLA-B*44:42",
                "HLA-B*44:43",
                "HLA-B*44:44",
                "HLA-B*44:45",
                "HLA-B*44:46",
                "HLA-B*44:47",
                "HLA-B*44:48",
                "HLA-B*44:49",
                "HLA-B*44:50",
                "HLA-B*44:51",
                "HLA-B*44:53",
                "HLA-B*44:54",
                "HLA-B*44:55",
                "HLA-B*44:57",
                "HLA-B*44:59",
                "HLA-B*44:60",
                "HLA-B*44:62",
                "HLA-B*44:63",
                "HLA-B*44:64",
                "HLA-B*44:65",
                "HLA-B*44:66",
                "HLA-B*44:67",
                "HLA-B*44:68",
                "HLA-B*44:69",
                "HLA-B*44:70",
                "HLA-B*44:71",
                "HLA-B*44:72",
                "HLA-B*44:73",
                "HLA-B*44:74",
                "HLA-B*44:75",
                "HLA-B*44:76",
                "HLA-B*44:77",
                "HLA-B*44:78",
                "HLA-B*44:79",
                "HLA-B*44:80",
                "HLA-B*44:81",
                "HLA-B*44:82",
                "HLA-B*44:83",
                "HLA-B*44:84",
                "HLA-B*44:85",
                "HLA-B*44:86",
                "HLA-B*44:87",
                "HLA-B*44:88",
                "HLA-B*44:89",
                "HLA-B*44:90",
                "HLA-B*44:91",
                "HLA-B*44:92",
                "HLA-B*44:93",
                "HLA-B*44:94",
                "HLA-B*44:95",
                "HLA-B*44:96",
                "HLA-B*44:97",
                "HLA-B*44:98",
                "HLA-B*44:99",
                "HLA-B*45:01",
                "HLA-B*45:02",
                "HLA-B*45:03",
                "HLA-B*45:04",
                "HLA-B*45:05",
                "HLA-B*45:06",
                "HLA-B*45:07",
                "HLA-B*45:08",
                "HLA-B*45:09",
                "HLA-B*45:10",
                "HLA-B*45:11",
                "HLA-B*45:12",
                "HLA-B*46:01",
                "HLA-B*46:02",
                "HLA-B*46:03",
                "HLA-B*46:04",
                "HLA-B*46:05",
                "HLA-B*46:06",
                "HLA-B*46:08",
                "HLA-B*46:09",
                "HLA-B*46:10",
                "HLA-B*46:11",
                "HLA-B*46:12",
                "HLA-B*46:13",
                "HLA-B*46:14",
                "HLA-B*46:16",
                "HLA-B*46:17",
                "HLA-B*46:18",
                "HLA-B*46:19",
                "HLA-B*46:20",
                "HLA-B*46:21",
                "HLA-B*46:22",
                "HLA-B*46:23",
                "HLA-B*46:24",
                "HLA-B*47:01",
                "HLA-B*47:02",
                "HLA-B*47:03",
                "HLA-B*47:04",
                "HLA-B*47:05",
                "HLA-B*47:06",
                "HLA-B*47:07",
                "HLA-B*48:01",
                "HLA-B*48:02",
                "HLA-B*48:03",
                "HLA-B*48:04",
                "HLA-B*48:05",
                "HLA-B*48:06",
                "HLA-B*48:07",
                "HLA-B*48:08",
                "HLA-B*48:09",
                "HLA-B*48:10",
                "HLA-B*48:11",
                "HLA-B*48:12",
                "HLA-B*48:13",
                "HLA-B*48:14",
                "HLA-B*48:15",
                "HLA-B*48:16",
                "HLA-B*48:17",
                "HLA-B*48:18",
                "HLA-B*48:19",
                "HLA-B*48:20",
                "HLA-B*48:21",
                "HLA-B*48:22",
                "HLA-B*48:23",
                "HLA-B*49:01",
                "HLA-B*49:02",
                "HLA-B*49:03",
                "HLA-B*49:04",
                "HLA-B*49:05",
                "HLA-B*49:06",
                "HLA-B*49:07",
                "HLA-B*49:08",
                "HLA-B*49:09",
                "HLA-B*49:10",
                "HLA-B*50:01",
                "HLA-B*50:02",
                "HLA-B*50:04",
                "HLA-B*50:05",
                "HLA-B*50:06",
                "HLA-B*50:07",
                "HLA-B*50:08",
                "HLA-B*50:09",
                "HLA-B*51:01",
                "HLA-B*51:02",
                "HLA-B*51:03",
                "HLA-B*51:04",
                "HLA-B*51:05",
                "HLA-B*51:06",
                "HLA-B*51:07",
                "HLA-B*51:08",
                "HLA-B*51:09",
                "HLA-B*51:12",
                "HLA-B*51:13",
                "HLA-B*51:14",
                "HLA-B*51:15",
                "HLA-B*51:16",
                "HLA-B*51:17",
                "HLA-B*51:18",
                "HLA-B*51:19",
                "HLA-B*51:20",
                "HLA-B*51:21",
                "HLA-B*51:22",
                "HLA-B*51:23",
                "HLA-B*51:24",
                "HLA-B*51:26",
                "HLA-B*51:28",
                "HLA-B*51:29",
                "HLA-B*51:30",
                "HLA-B*51:31",
                "HLA-B*51:32",
                "HLA-B*51:33",
                "HLA-B*51:34",
                "HLA-B*51:35",
                "HLA-B*51:36",
                "HLA-B*51:37",
                "HLA-B*51:38",
                "HLA-B*51:39",
                "HLA-B*51:40",
                "HLA-B*51:42",
                "HLA-B*51:43",
                "HLA-B*51:45",
                "HLA-B*51:46",
                "HLA-B*51:48",
                "HLA-B*51:49",
                "HLA-B*51:50",
                "HLA-B*51:51",
                "HLA-B*51:52",
                "HLA-B*51:53",
                "HLA-B*51:54",
                "HLA-B*51:55",
                "HLA-B*51:56",
                "HLA-B*51:57",
                "HLA-B*51:58",
                "HLA-B*51:59",
                "HLA-B*51:60",
                "HLA-B*51:61",
                "HLA-B*51:62",
                "HLA-B*51:63",
                "HLA-B*51:64",
                "HLA-B*51:65",
                "HLA-B*51:66",
                "HLA-B*51:67",
                "HLA-B*51:68",
                "HLA-B*51:69",
                "HLA-B*51:70",
                "HLA-B*51:71",
                "HLA-B*51:72",
                "HLA-B*51:73",
                "HLA-B*51:74",
                "HLA-B*51:75",
                "HLA-B*51:76",
                "HLA-B*51:77",
                "HLA-B*51:78",
                "HLA-B*51:79",
                "HLA-B*51:80",
                "HLA-B*51:81",
                "HLA-B*51:82",
                "HLA-B*51:83",
                "HLA-B*51:84",
                "HLA-B*51:85",
                "HLA-B*51:86",
                "HLA-B*51:87",
                "HLA-B*51:88",
                "HLA-B*51:89",
                "HLA-B*51:90",
                "HLA-B*51:91",
                "HLA-B*51:92",
                "HLA-B*51:93",
                "HLA-B*51:94",
                "HLA-B*51:95",
                "HLA-B*51:96",
                "HLA-B*52:01",
                "HLA-B*52:02",
                "HLA-B*52:03",
                "HLA-B*52:04",
                "HLA-B*52:05",
                "HLA-B*52:06",
                "HLA-B*52:07",
                "HLA-B*52:08",
                "HLA-B*52:09",
                "HLA-B*52:10",
                "HLA-B*52:11",
                "HLA-B*52:12",
                "HLA-B*52:13",
                "HLA-B*52:14",
                "HLA-B*52:15",
                "HLA-B*52:16",
                "HLA-B*52:17",
                "HLA-B*52:18",
                "HLA-B*52:19",
                "HLA-B*52:20",
                "HLA-B*52:21",
                "HLA-B*53:01",
                "HLA-B*53:02",
                "HLA-B*53:03",
                "HLA-B*53:04",
                "HLA-B*53:05",
                "HLA-B*53:06",
                "HLA-B*53:07",
                "HLA-B*53:08",
                "HLA-B*53:09",
                "HLA-B*53:10",
                "HLA-B*53:11",
                "HLA-B*53:12",
                "HLA-B*53:13",
                "HLA-B*53:14",
                "HLA-B*53:15",
                "HLA-B*53:16",
                "HLA-B*53:17",
                "HLA-B*53:18",
                "HLA-B*53:19",
                "HLA-B*53:20",
                "HLA-B*53:21",
                "HLA-B*53:22",
                "HLA-B*53:23",
                "HLA-B*54:01",
                "HLA-B*54:02",
                "HLA-B*54:03",
                "HLA-B*54:04",
                "HLA-B*54:06",
                "HLA-B*54:07",
                "HLA-B*54:09",
                "HLA-B*54:10",
                "HLA-B*54:11",
                "HLA-B*54:12",
                "HLA-B*54:13",
                "HLA-B*54:14",
                "HLA-B*54:15",
                "HLA-B*54:16",
                "HLA-B*54:17",
                "HLA-B*54:18",
                "HLA-B*54:19",
                "HLA-B*54:20",
                "HLA-B*54:21",
                "HLA-B*54:22",
                "HLA-B*54:23",
                "HLA-B*55:01",
                "HLA-B*55:02",
                "HLA-B*55:03",
                "HLA-B*55:04",
                "HLA-B*55:05",
                "HLA-B*55:07",
                "HLA-B*55:08",
                "HLA-B*55:09",
                "HLA-B*55:10",
                "HLA-B*55:11",
                "HLA-B*55:12",
                "HLA-B*55:13",
                "HLA-B*55:14",
                "HLA-B*55:15",
                "HLA-B*55:16",
                "HLA-B*55:17",
                "HLA-B*55:18",
                "HLA-B*55:19",
                "HLA-B*55:20",
                "HLA-B*55:21",
                "HLA-B*55:22",
                "HLA-B*55:23",
                "HLA-B*55:24",
                "HLA-B*55:25",
                "HLA-B*55:26",
                "HLA-B*55:27",
                "HLA-B*55:28",
                "HLA-B*55:29",
                "HLA-B*55:30",
                "HLA-B*55:31",
                "HLA-B*55:32",
                "HLA-B*55:33",
                "HLA-B*55:34",
                "HLA-B*55:35",
                "HLA-B*55:36",
                "HLA-B*55:37",
                "HLA-B*55:38",
                "HLA-B*55:39",
                "HLA-B*55:40",
                "HLA-B*55:41",
                "HLA-B*55:42",
                "HLA-B*55:43",
                "HLA-B*56:01",
                "HLA-B*56:02",
                "HLA-B*56:03",
                "HLA-B*56:04",
                "HLA-B*56:05",
                "HLA-B*56:06",
                "HLA-B*56:07",
                "HLA-B*56:08",
                "HLA-B*56:09",
                "HLA-B*56:10",
                "HLA-B*56:11",
                "HLA-B*56:12",
                "HLA-B*56:13",
                "HLA-B*56:14",
                "HLA-B*56:15",
                "HLA-B*56:16",
                "HLA-B*56:17",
                "HLA-B*56:18",
                "HLA-B*56:20",
                "HLA-B*56:21",
                "HLA-B*56:22",
                "HLA-B*56:23",
                "HLA-B*56:24",
                "HLA-B*56:25",
                "HLA-B*56:26",
                "HLA-B*56:27",
                "HLA-B*56:29",
                "HLA-B*57:01",
                "HLA-B*57:02",
                "HLA-B*57:03",
                "HLA-B*57:04",
                "HLA-B*57:05",
                "HLA-B*57:06",
                "HLA-B*57:07",
                "HLA-B*57:08",
                "HLA-B*57:09",
                "HLA-B*57:10",
                "HLA-B*57:11",
                "HLA-B*57:12",
                "HLA-B*57:13",
                "HLA-B*57:14",
                "HLA-B*57:15",
                "HLA-B*57:16",
                "HLA-B*57:17",
                "HLA-B*57:18",
                "HLA-B*57:19",
                "HLA-B*57:20",
                "HLA-B*57:21",
                "HLA-B*57:22",
                "HLA-B*57:23",
                "HLA-B*57:24",
                "HLA-B*57:25",
                "HLA-B*57:26",
                "HLA-B*57:27",
                "HLA-B*57:29",
                "HLA-B*57:30",
                "HLA-B*57:31",
                "HLA-B*57:32",
                "HLA-B*58:01",
                "HLA-B*58:02",
                "HLA-B*58:04",
                "HLA-B*58:05",
                "HLA-B*58:06",
                "HLA-B*58:07",
                "HLA-B*58:08",
                "HLA-B*58:09",
                "HLA-B*58:11",
                "HLA-B*58:12",
                "HLA-B*58:13",
                "HLA-B*58:14",
                "HLA-B*58:15",
                "HLA-B*58:16",
                "HLA-B*58:18",
                "HLA-B*58:19",
                "HLA-B*58:20",
                "HLA-B*58:21",
                "HLA-B*58:22",
                "HLA-B*58:23",
                "HLA-B*58:24",
                "HLA-B*58:25",
                "HLA-B*58:26",
                "HLA-B*58:27",
                "HLA-B*58:28",
                "HLA-B*58:29",
                "HLA-B*58:30",
                "HLA-B*59:01",
                "HLA-B*59:02",
                "HLA-B*59:03",
                "HLA-B*59:04",
                "HLA-B*59:05",
                "HLA-B*67:01",
                "HLA-B*67:02",
                "HLA-B*73:01",
                "HLA-B*73:02",
                "HLA-B*78:01",
                "HLA-B*78:02",
                "HLA-B*78:03",
                "HLA-B*78:04",
                "HLA-B*78:05",
                "HLA-B*78:06",
                "HLA-B*78:07",
                "HLA-B*81:01",
                "HLA-B*81:02",
                "HLA-B*81:03",
                "HLA-B*81:05",
                "HLA-B*82:01",
                "HLA-B*82:02",
                "HLA-B*82:03",
                "HLA-B*83:01",
                "HLA-C*01:02",
                "HLA-C*01:03",
                "HLA-C*01:04",
                "HLA-C*01:05",
                "HLA-C*01:06",
                "HLA-C*01:07",
                "HLA-C*01:08",
                "HLA-C*01:09",
                "HLA-C*01:10",
                "HLA-C*01:11",
                "HLA-C*01:12",
                "HLA-C*01:13",
                "HLA-C*01:14",
                "HLA-C*01:15",
                "HLA-C*01:16",
                "HLA-C*01:17",
                "HLA-C*01:18",
                "HLA-C*01:19",
                "HLA-C*01:20",
                "HLA-C*01:21",
                "HLA-C*01:22",
                "HLA-C*01:23",
                "HLA-C*01:24",
                "HLA-C*01:25",
                "HLA-C*01:26",
                "HLA-C*01:27",
                "HLA-C*01:28",
                "HLA-C*01:29",
                "HLA-C*01:30",
                "HLA-C*01:31",
                "HLA-C*01:32",
                "HLA-C*01:33",
                "HLA-C*01:34",
                "HLA-C*01:35",
                "HLA-C*01:36",
                "HLA-C*01:38",
                "HLA-C*01:39",
                "HLA-C*01:40",
                "HLA-C*02:02",
                "HLA-C*02:03",
                "HLA-C*02:04",
                "HLA-C*02:05",
                "HLA-C*02:06",
                "HLA-C*02:07",
                "HLA-C*02:08",
                "HLA-C*02:09",
                "HLA-C*02:10",
                "HLA-C*02:11",
                "HLA-C*02:12",
                "HLA-C*02:13",
                "HLA-C*02:14",
                "HLA-C*02:15",
                "HLA-C*02:16",
                "HLA-C*02:17",
                "HLA-C*02:18",
                "HLA-C*02:19",
                "HLA-C*02:20",
                "HLA-C*02:21",
                "HLA-C*02:22",
                "HLA-C*02:23",
                "HLA-C*02:24",
                "HLA-C*02:26",
                "HLA-C*02:27",
                "HLA-C*02:28",
                "HLA-C*02:29",
                "HLA-C*02:30",
                "HLA-C*02:31",
                "HLA-C*02:32",
                "HLA-C*02:33",
                "HLA-C*02:34",
                "HLA-C*02:35",
                "HLA-C*02:36",
                "HLA-C*02:37",
                "HLA-C*02:39",
                "HLA-C*02:40",
                "HLA-C*03:01",
                "HLA-C*03:02",
                "HLA-C*03:03",
                "HLA-C*03:04",
                "HLA-C*03:05",
                "HLA-C*03:06",
                "HLA-C*03:07",
                "HLA-C*03:08",
                "HLA-C*03:09",
                "HLA-C*03:10",
                "HLA-C*03:11",
                "HLA-C*03:12",
                "HLA-C*03:13",
                "HLA-C*03:14",
                "HLA-C*03:15",
                "HLA-C*03:16",
                "HLA-C*03:17",
                "HLA-C*03:18",
                "HLA-C*03:19",
                "HLA-C*03:21",
                "HLA-C*03:23",
                "HLA-C*03:24",
                "HLA-C*03:25",
                "HLA-C*03:26",
                "HLA-C*03:27",
                "HLA-C*03:28",
                "HLA-C*03:29",
                "HLA-C*03:30",
                "HLA-C*03:31",
                "HLA-C*03:32",
                "HLA-C*03:33",
                "HLA-C*03:34",
                "HLA-C*03:35",
                "HLA-C*03:36",
                "HLA-C*03:37",
                "HLA-C*03:38",
                "HLA-C*03:39",
                "HLA-C*03:40",
                "HLA-C*03:41",
                "HLA-C*03:42",
                "HLA-C*03:43",
                "HLA-C*03:44",
                "HLA-C*03:45",
                "HLA-C*03:46",
                "HLA-C*03:47",
                "HLA-C*03:48",
                "HLA-C*03:49",
                "HLA-C*03:50",
                "HLA-C*03:51",
                "HLA-C*03:52",
                "HLA-C*03:53",
                "HLA-C*03:54",
                "HLA-C*03:55",
                "HLA-C*03:56",
                "HLA-C*03:57",
                "HLA-C*03:58",
                "HLA-C*03:59",
                "HLA-C*03:60",
                "HLA-C*03:61",
                "HLA-C*03:62",
                "HLA-C*03:63",
                "HLA-C*03:64",
                "HLA-C*03:65",
                "HLA-C*03:66",
                "HLA-C*03:67",
                "HLA-C*03:68",
                "HLA-C*03:69",
                "HLA-C*03:70",
                "HLA-C*03:71",
                "HLA-C*03:72",
                "HLA-C*03:73",
                "HLA-C*03:74",
                "HLA-C*03:75",
                "HLA-C*03:76",
                "HLA-C*03:77",
                "HLA-C*03:78",
                "HLA-C*03:79",
                "HLA-C*03:80",
                "HLA-C*03:81",
                "HLA-C*03:82",
                "HLA-C*03:83",
                "HLA-C*03:84",
                "HLA-C*03:85",
                "HLA-C*03:86",
                "HLA-C*03:87",
                "HLA-C*03:88",
                "HLA-C*03:89",
                "HLA-C*03:90",
                "HLA-C*03:91",
                "HLA-C*03:92",
                "HLA-C*03:93",
                "HLA-C*03:94",
                "HLA-C*04:01",
                "HLA-C*04:03",
                "HLA-C*04:04",
                "HLA-C*04:05",
                "HLA-C*04:06",
                "HLA-C*04:07",
                "HLA-C*04:08",
                "HLA-C*04:10",
                "HLA-C*04:11",
                "HLA-C*04:12",
                "HLA-C*04:13",
                "HLA-C*04:14",
                "HLA-C*04:15",
                "HLA-C*04:16",
                "HLA-C*04:17",
                "HLA-C*04:18",
                "HLA-C*04:19",
                "HLA-C*04:20",
                "HLA-C*04:23",
                "HLA-C*04:24",
                "HLA-C*04:25",
                "HLA-C*04:26",
                "HLA-C*04:27",
                "HLA-C*04:28",
                "HLA-C*04:29",
                "HLA-C*04:30",
                "HLA-C*04:31",
                "HLA-C*04:32",
                "HLA-C*04:33",
                "HLA-C*04:34",
                "HLA-C*04:35",
                "HLA-C*04:36",
                "HLA-C*04:37",
                "HLA-C*04:38",
                "HLA-C*04:39",
                "HLA-C*04:40",
                "HLA-C*04:41",
                "HLA-C*04:42",
                "HLA-C*04:43",
                "HLA-C*04:44",
                "HLA-C*04:45",
                "HLA-C*04:46",
                "HLA-C*04:47",
                "HLA-C*04:48",
                "HLA-C*04:49",
                "HLA-C*04:50",
                "HLA-C*04:51",
                "HLA-C*04:52",
                "HLA-C*04:53",
                "HLA-C*04:54",
                "HLA-C*04:55",
                "HLA-C*04:56",
                "HLA-C*04:57",
                "HLA-C*04:58",
                "HLA-C*04:60",
                "HLA-C*04:61",
                "HLA-C*04:62",
                "HLA-C*04:63",
                "HLA-C*04:64",
                "HLA-C*04:65",
                "HLA-C*04:66",
                "HLA-C*04:67",
                "HLA-C*04:68",
                "HLA-C*04:69",
                "HLA-C*04:70",
                "HLA-C*05:01",
                "HLA-C*05:03",
                "HLA-C*05:04",
                "HLA-C*05:05",
                "HLA-C*05:06",
                "HLA-C*05:08",
                "HLA-C*05:09",
                "HLA-C*05:10",
                "HLA-C*05:11",
                "HLA-C*05:12",
                "HLA-C*05:13",
                "HLA-C*05:14",
                "HLA-C*05:15",
                "HLA-C*05:16",
                "HLA-C*05:17",
                "HLA-C*05:18",
                "HLA-C*05:19",
                "HLA-C*05:20",
                "HLA-C*05:21",
                "HLA-C*05:22",
                "HLA-C*05:23",
                "HLA-C*05:24",
                "HLA-C*05:25",
                "HLA-C*05:26",
                "HLA-C*05:27",
                "HLA-C*05:28",
                "HLA-C*05:29",
                "HLA-C*05:30",
                "HLA-C*05:31",
                "HLA-C*05:32",
                "HLA-C*05:33",
                "HLA-C*05:34",
                "HLA-C*05:35",
                "HLA-C*05:36",
                "HLA-C*05:37",
                "HLA-C*05:38",
                "HLA-C*05:39",
                "HLA-C*05:40",
                "HLA-C*05:41",
                "HLA-C*05:42",
                "HLA-C*05:43",
                "HLA-C*05:44",
                "HLA-C*05:45",
                "HLA-C*06:02",
                "HLA-C*06:03",
                "HLA-C*06:04",
                "HLA-C*06:05",
                "HLA-C*06:06",
                "HLA-C*06:07",
                "HLA-C*06:08",
                "HLA-C*06:09",
                "HLA-C*06:10",
                "HLA-C*06:11",
                "HLA-C*06:12",
                "HLA-C*06:13",
                "HLA-C*06:14",
                "HLA-C*06:15",
                "HLA-C*06:17",
                "HLA-C*06:18",
                "HLA-C*06:19",
                "HLA-C*06:20",
                "HLA-C*06:21",
                "HLA-C*06:22",
                "HLA-C*06:23",
                "HLA-C*06:24",
                "HLA-C*06:25",
                "HLA-C*06:26",
                "HLA-C*06:27",
                "HLA-C*06:28",
                "HLA-C*06:29",
                "HLA-C*06:30",
                "HLA-C*06:31",
                "HLA-C*06:32",
                "HLA-C*06:33",
                "HLA-C*06:34",
                "HLA-C*06:35",
                "HLA-C*06:36",
                "HLA-C*06:37",
                "HLA-C*06:38",
                "HLA-C*06:39",
                "HLA-C*06:40",
                "HLA-C*06:41",
                "HLA-C*06:42",
                "HLA-C*06:43",
                "HLA-C*06:44",
                "HLA-C*06:45",
                "HLA-C*07:01",
                "HLA-C*07:02",
                "HLA-C*07:03",
                "HLA-C*07:04",
                "HLA-C*07:05",
                "HLA-C*07:06",
                "HLA-C*07:07",
                "HLA-C*07:08",
                "HLA-C*07:09",
                "HLA-C*07:10",
                "HLA-C*07:100",
                "HLA-C*07:101",
                "HLA-C*07:102",
                "HLA-C*07:103",
                "HLA-C*07:105",
                "HLA-C*07:106",
                "HLA-C*07:107",
                "HLA-C*07:108",
                "HLA-C*07:109",
                "HLA-C*07:11",
                "HLA-C*07:110",
                "HLA-C*07:111",
                "HLA-C*07:112",
                "HLA-C*07:113",
                "HLA-C*07:114",
                "HLA-C*07:115",
                "HLA-C*07:116",
                "HLA-C*07:117",
                "HLA-C*07:118",
                "HLA-C*07:119",
                "HLA-C*07:12",
                "HLA-C*07:120",
                "HLA-C*07:122",
                "HLA-C*07:123",
                "HLA-C*07:124",
                "HLA-C*07:125",
                "HLA-C*07:126",
                "HLA-C*07:127",
                "HLA-C*07:128",
                "HLA-C*07:129",
                "HLA-C*07:13",
                "HLA-C*07:130",
                "HLA-C*07:131",
                "HLA-C*07:132",
                "HLA-C*07:133",
                "HLA-C*07:134",
                "HLA-C*07:135",
                "HLA-C*07:136",
                "HLA-C*07:137",
                "HLA-C*07:138",
                "HLA-C*07:139",
                "HLA-C*07:14",
                "HLA-C*07:140",
                "HLA-C*07:141",
                "HLA-C*07:142",
                "HLA-C*07:143",
                "HLA-C*07:144",
                "HLA-C*07:145",
                "HLA-C*07:146",
                "HLA-C*07:147",
                "HLA-C*07:148",
                "HLA-C*07:149",
                "HLA-C*07:15",
                "HLA-C*07:16",
                "HLA-C*07:17",
                "HLA-C*07:18",
                "HLA-C*07:19",
                "HLA-C*07:20",
                "HLA-C*07:21",
                "HLA-C*07:22",
                "HLA-C*07:23",
                "HLA-C*07:24",
                "HLA-C*07:25",
                "HLA-C*07:26",
                "HLA-C*07:27",
                "HLA-C*07:28",
                "HLA-C*07:29",
                "HLA-C*07:30",
                "HLA-C*07:31",
                "HLA-C*07:35",
                "HLA-C*07:36",
                "HLA-C*07:37",
                "HLA-C*07:38",
                "HLA-C*07:39",
                "HLA-C*07:40",
                "HLA-C*07:41",
                "HLA-C*07:42",
                "HLA-C*07:43",
                "HLA-C*07:44",
                "HLA-C*07:45",
                "HLA-C*07:46",
                "HLA-C*07:47",
                "HLA-C*07:48",
                "HLA-C*07:49",
                "HLA-C*07:50",
                "HLA-C*07:51",
                "HLA-C*07:52",
                "HLA-C*07:53",
                "HLA-C*07:54",
                "HLA-C*07:56",
                "HLA-C*07:57",
                "HLA-C*07:58",
                "HLA-C*07:59",
                "HLA-C*07:60",
                "HLA-C*07:62",
                "HLA-C*07:63",
                "HLA-C*07:64",
                "HLA-C*07:65",
                "HLA-C*07:66",
                "HLA-C*07:67",
                "HLA-C*07:68",
                "HLA-C*07:69",
                "HLA-C*07:70",
                "HLA-C*07:71",
                "HLA-C*07:72",
                "HLA-C*07:73",
                "HLA-C*07:74",
                "HLA-C*07:75",
                "HLA-C*07:76",
                "HLA-C*07:77",
                "HLA-C*07:78",
                "HLA-C*07:79",
                "HLA-C*07:80",
                "HLA-C*07:81",
                "HLA-C*07:82",
                "HLA-C*07:83",
                "HLA-C*07:84",
                "HLA-C*07:85",
                "HLA-C*07:86",
                "HLA-C*07:87",
                "HLA-C*07:88",
                "HLA-C*07:89",
                "HLA-C*07:90",
                "HLA-C*07:91",
                "HLA-C*07:92",
                "HLA-C*07:93",
                "HLA-C*07:94",
                "HLA-C*07:95",
                "HLA-C*07:96",
                "HLA-C*07:97",
                "HLA-C*07:99",
                "HLA-C*08:01",
                "HLA-C*08:02",
                "HLA-C*08:03",
                "HLA-C*08:04",
                "HLA-C*08:05",
                "HLA-C*08:06",
                "HLA-C*08:07",
                "HLA-C*08:08",
                "HLA-C*08:09",
                "HLA-C*08:10",
                "HLA-C*08:11",
                "HLA-C*08:12",
                "HLA-C*08:13",
                "HLA-C*08:14",
                "HLA-C*08:15",
                "HLA-C*08:16",
                "HLA-C*08:17",
                "HLA-C*08:18",
                "HLA-C*08:19",
                "HLA-C*08:20",
                "HLA-C*08:21",
                "HLA-C*08:22",
                "HLA-C*08:23",
                "HLA-C*08:24",
                "HLA-C*08:25",
                "HLA-C*08:27",
                "HLA-C*08:28",
                "HLA-C*08:29",
                "HLA-C*08:30",
                "HLA-C*08:31",
                "HLA-C*08:32",
                "HLA-C*08:33",
                "HLA-C*08:34",
                "HLA-C*08:35",
                "HLA-C*12:02",
                "HLA-C*12:03",
                "HLA-C*12:04",
                "HLA-C*12:05",
                "HLA-C*12:06",
                "HLA-C*12:07",
                "HLA-C*12:08",
                "HLA-C*12:09",
                "HLA-C*12:10",
                "HLA-C*12:11",
                "HLA-C*12:12",
                "HLA-C*12:13",
                "HLA-C*12:14",
                "HLA-C*12:15",
                "HLA-C*12:16",
                "HLA-C*12:17",
                "HLA-C*12:18",
                "HLA-C*12:19",
                "HLA-C*12:20",
                "HLA-C*12:21",
                "HLA-C*12:22",
                "HLA-C*12:23",
                "HLA-C*12:24",
                "HLA-C*12:25",
                "HLA-C*12:26",
                "HLA-C*12:27",
                "HLA-C*12:28",
                "HLA-C*12:29",
                "HLA-C*12:30",
                "HLA-C*12:31",
                "HLA-C*12:32",
                "HLA-C*12:33",
                "HLA-C*12:34",
                "HLA-C*12:35",
                "HLA-C*12:36",
                "HLA-C*12:37",
                "HLA-C*12:38",
                "HLA-C*12:40",
                "HLA-C*12:41",
                "HLA-C*12:43",
                "HLA-C*12:44",
                "HLA-C*14:02",
                "HLA-C*14:03",
                "HLA-C*14:04",
                "HLA-C*14:05",
                "HLA-C*14:06",
                "HLA-C*14:08",
                "HLA-C*14:09",
                "HLA-C*14:10",
                "HLA-C*14:11",
                "HLA-C*14:12",
                "HLA-C*14:13",
                "HLA-C*14:14",
                "HLA-C*14:15",
                "HLA-C*14:16",
                "HLA-C*14:17",
                "HLA-C*14:18",
                "HLA-C*14:19",
                "HLA-C*14:20",
                "HLA-C*15:02",
                "HLA-C*15:03",
                "HLA-C*15:04",
                "HLA-C*15:05",
                "HLA-C*15:06",
                "HLA-C*15:07",
                "HLA-C*15:08",
                "HLA-C*15:09",
                "HLA-C*15:10",
                "HLA-C*15:11",
                "HLA-C*15:12",
                "HLA-C*15:13",
                "HLA-C*15:15",
                "HLA-C*15:16",
                "HLA-C*15:17",
                "HLA-C*15:18",
                "HLA-C*15:19",
                "HLA-C*15:20",
                "HLA-C*15:21",
                "HLA-C*15:22",
                "HLA-C*15:23",
                "HLA-C*15:24",
                "HLA-C*15:25",
                "HLA-C*15:26",
                "HLA-C*15:27",
                "HLA-C*15:28",
                "HLA-C*15:29",
                "HLA-C*15:30",
                "HLA-C*15:31",
                "HLA-C*15:33",
                "HLA-C*15:34",
                "HLA-C*15:35",
                "HLA-C*16:01",
                "HLA-C*16:02",
                "HLA-C*16:04",
                "HLA-C*16:06",
                "HLA-C*16:07",
                "HLA-C*16:08",
                "HLA-C*16:09",
                "HLA-C*16:10",
                "HLA-C*16:11",
                "HLA-C*16:12",
                "HLA-C*16:13",
                "HLA-C*16:14",
                "HLA-C*16:15",
                "HLA-C*16:17",
                "HLA-C*16:18",
                "HLA-C*16:19",
                "HLA-C*16:20",
                "HLA-C*16:21",
                "HLA-C*16:22",
                "HLA-C*16:23",
                "HLA-C*16:24",
                "HLA-C*16:25",
                "HLA-C*16:26",
                "HLA-C*17:01",
                "HLA-C*17:02",
                "HLA-C*17:03",
                "HLA-C*17:04",
                "HLA-C*17:05",
                "HLA-C*17:06",
                "HLA-C*17:07",
                "HLA-C*18:01",
                "HLA-C*18:02",
                "HLA-C*18:03",
                "HLA-E*01:01",
                "HLA-G*01:01",
                "HLA-G*01:02",
                "HLA-G*01:03",
                "HLA-G*01:04",
                "HLA-G*01:06",
                "HLA-G*01:07",
                "HLA-G*01:08",
                "HLA-G*01:09",
                "H-2-Db",
                "H-2-Dd",
                "H-2-Kb",
                "H-2-Kd",
                "H-2-Kk",
                "H-2-Ld"
    ]

    if allele in allele_list:
        return True
    else:
        return False


main()


#def determine_mt_pep_seq(mt_seq, 



            

                # for t in transcripts:

                    # csq = resolve_consequence(t['Consequence'])
                    # if csq is None:
                        # continue

                    # elif csq == "frameshift":
                        # if t["NMD"] == "":
                            # continue
                        # #elif tranascript['NMD'] == 'NMD_escaping_variant':
                    # else:
                        # amino_acid_change_position = (
                            # t["Protein_position"] + '_' + t["Amino_acids"]
                        # )

                    # #print(alt)

                    # # determine amino acid sequences (MT,WT)
                    # wt_amino_acid_seq = pep_mod[t_name]
                    # mt_amino_acid_seq = determine_mt_peptide(
                            # wt_amino_acid_seq, 
                            # t['Protein_position'],
                            # t['Amino_acids'])
                    # mt_pep_seq = determine_subpeptide(
                            # mt_amino_acid_seq, 
                            # t['Protein_position'],
                            # t["Amino_acids"])


                    # gene_name = t["SYMBOL"]
                    # #index = f"{gene_name}_{t_name}_{transcript_count[t_name]}.{csq}.{amino_acid_change_position}"
                    # ensembl_gene_id = t["Gene"]
                    # output_row = {
                        # "chrom": entry.CHROM,
                        # "start": entry.affected_start,
                        # "stop": entry.affected_end,
                        # "source": entry.INFO["SRC"],
                        # "group": entry.INFO["GRP"],
                        # "reference": entry.REF,
                        # "variant": alt.value,
                        # "gene_name": gene_name,
                        # "transcript_name": t_name,
                        # "amino_acid_change": t["Amino_acids"],
                        # "ensembl_gene_id": ensembl_gene_id,
                        # "wildtype_amino_acid_sequence": pep_mod[t_name],
                        # "mutant_amino_acid_sequence": mt_amino_acid_seq,
                        # "mutant_peptide_sequence": mt_pep_seq,
                        # "variant_type": csq,
                        # "protein_position": t["Protein_position"]
                    # }





# #                    iedb_call(





                    # for x in keys:
                        # outputfile.write(str(output_row[x]) + '\t')
                    # outputfile.write('\n')



#                if alt.type == 'DEL':


#                    bam_readcount_position = start + 1
                # if alt.type == 'DEL':
                    # (simpl_ref, simpl_alt) = simplify_indel_allele(ref, alt)
                    # ref_base = ref[1:2]
                    # var_base = "-" + simpl_ref
                    # variant_type = "indels"

                # if alt.type == 'INS':
                    # bam_readcount_position = start
                    # (simpl_ref, simpl_alt) = simplify_indel_allele(ref, alt)
                    # ref_base = ref
                    # var_base = "+" + simpl_alt
                    # variant_type = "indels"


                # if alt.type == 'SNV':
                    # print("pos: " + str(entry.POS))
                    # print("ref: " + str(entry.REF))
                    # # print(alt.type)
                    # # print(entry)
                    # # print(genome[entry.CHROM][entry.POS-1])
                # elif alt.type == 'DEL':
                    # print("pos: " + str(entry.POS))
                    # print("ref: " + str(entry.REF))
                    # print(alt.type)
                    # print("alt: " + str(alt))
                    # print(genome[entry.CHROM][entry.POS-2:entry.POS+2])
                    # print(entry.affected_start)
                    # print(entry.affected_end)
            
                    # print(alleles_dict)

                # # elif alt.type == 'INS':
                    # # print("pos: " + str(entry.POS))
                    # # print("ref: " + str(entry.REF))
                    # # print(alt.type)
                    # # print("alt: " + str(alt))
                    # # print(genome[entry.CHROM][entry.POS-2:entry.POS+2])
                    # # print(entry.affected_start)
                    # # print(entry.affected_end)
            





                    # amino_acid_change_position = transcript["Protein_position"]
                    # amino_acid_change_position = (
                        # transcript["Protein_position"] + transcript["Amino_acids"]
                    # )



                # if csq is None:
                    # continue
                # else:
                    # print(csq)



                # elif csq == 'frameshift':
                    # print("frameshift")



                # if name != '':
                    # print(name)
                        



                                        

                        #output_row["mutant_peptide_sequence"] = field["DownstreamProtein"]













                    # print(csq)
                    # print(f"Protein position {field['Protein_position']}")
                    # print(f"Amino Acids {field['Amino_acids']}")
                    # print(f"Downstream {field['DownstreamProtein']}")
                    # print(field["Codons"])
                    # print(field["NMD"])
                    # print(f"WildtypeProtein: {field['WildtypeProtein']}")
                    
                    # # determine amino acid sequences (MT,WT)
                    # wt_aa_seq = pep_mod[name]  # determines the WT sequence
                    # print(wt_aa_seq)
                    # print(f"wt: {wt_aa_seq}")
                    # mt_aa_seq = determine_mt_aa_seq(
                            # csq,
                            # wt_aa_seq, 
                            # field['Protein_position'],
                            # field['Amino_acids'])
                    # print(f"mt: {mt_aa_seq}")
                            





                            # var_start_pos = get_variant_startpos(output_row['aa_change'])
                            # print("startpos")




                            # print(f"position: {pos}")

                            # # get subsequence
                            # (mt_start, wt_subseq, mt_subseq) = get_frameshift_subseq(
                                    # pos,
                                    # output_row["wt_aa_seq"],
                                    # epitope_len,
                                    # output_row
                            # )

                            # print(f"wt_subseq: {wt_subseq}")
                            # mt_subseq += output_row['downstream_aa_seq']
                            # print(f"mt_subseq: {mt_subseq}")


                        # elif (output_row["variant_type"] == "missense" 
                              # or output_row["variant_type"] == "inframe_ins" 
                              # or output_row["variant_type"] == "inframe_del"):
                            # wt_aa, mt_aa = output_row['aa_change'].split('/')

                            # # consider stop codons
                            # wt_aa, mt_aa, stop_codon = scan_stop_codons(wt_aa, mt_aa)

                            # if (output_row["variant_type"] == "missense" or 
                                # output_row["variant_type"] == "inframe_ins"):
                                # if "-" in output_row["protein_position"]:
                                    # pos = int(output_row["protein_position"].split('-',1)[0]) - 1
                                    # wt_aa_len = len(wt_aa)
                                # else:  # single variant
                                    # pos = int(output_row["protein_position"]) - 1
                                    # wt_aa_len = len(wt_aa)

                            # if output_row["variant_type"] == "inframe_del":
                                # pos = int(output_row["protein_position"].split('-',1)[0]) - 1
                                # wt_aa_len = len(wt_aa)

                            # if mt_aa == "-":
                                # mt_aa = ''

                            # mt_start, wt_subseq = get_wt_subseq(
                                    # pos,
                                    # output_row['wt_aa_seq'],
                                    # wt_aa_len,
                                    # epitope_len,
                                    # output_row)
                            # mt_end = mt_start + wt_aa_len
                            # if (wt_aa != '-' and wt_aa != wt_subseq[mt_start:mt_end]):
                                # continue

                            # if stop_codon:
                                # mt_subseq = wt_subseq[:mt_start] + mt_aa
                            # else:
                                # mt_subseq = wt_subseq[:mt_start] + mt_aa + wt_subseq[mt_end:]




                        # # final exclusions
                        # if '*' in wt_subseq or '*' in mt_subseq:
                            # continue
                        # if 'X' in wt_subseq or 'X' in mt_subseq:
                            # continue
                        # if mt_subseq in wt_subseq:
                            # continue

                        # if (len(wt_subseq) < epitope_len or
                            # len(mt_subseq) < epitope_len):
                            # continue



                        # calculate binding affinity
#                        with tempfile.TemporaryDirectory() as tmpdir:





#def generate_subsequence():

                        # end_pos_mt = start_pos_var + len(output_row['downstream_aa_seq'])-1  
                        # mt_start_pos = 

                        # endpos = startpos + len(output_row['downstream_aa_seq']) - 1
                        # wt_seq = output_row['wt_aa_seq']  # wildtype sequence
                        # mt_seq = wt_seq[:startpos] + output_row['downstream_aa_seq']

                        # wt_fasta.write(f">id{output_row['seqnum']}\n{wt_seq}\n")
                        # mt_fasta.write(f">{output_row['seqnum']}\n{mt_seq}\n")

                        # fill up with $ to have wt and mt of same length
                        # if len(mt_seq) > len(wt_seq):
                            # for i in range(len(wt_seq), len(mt_seq)+1):
                                # wt_seq += '$'
                        # if len(wt_seq) > len(mt_seq):
                            # for i in range(len(mt_seq), len(wt_seq)+1):
                                # mt_seq += '$'
                        # elif len(mt_seq) > len(wt_seq):
                            # for i in range(len(wt_seq), len(mt_seq)+1):
                                # wt_seq += '$'

                        # print(f"wt_seq: {wt_seq}")
                        # print(f"protein position: {output_row['protein_position']}")
                        # print(f"start position variant: {startpos}")
                        # print(f"end position variant: {endpos}")
                        # print(f"aa change: {output_row['aa_change']}")
                        # print(f"downstream: {output_row['downstream_aa_seq']}")
                        # print(f"mt_seq: {mt_seq}")

                        # print(f"epitope_len: {epitope_length}")
                        # wt_subseq, mt_subseq = determine_fs_subsequences(
                            # wt_seq,
                            # mt_seq,
                            # epitope_length,
                            # startpos,
                            # endpos
                        # )

                        # print(wt_subseq)
                        # print(mt_subseq)

#                    wt_fasta.write('>id'



                        






                        




                        # for key in subseqs.keys():
                            # output_row['wt_pep'] = subseqs[key][0]
                            # output_row['mt_pep'] = subseqs[key][1]

                            # out_itrm.append(output_row)
    # print(out_itrm)

    # calculate binding affinity
    #wt_binding_affinity = calc_peptide_binding(alleles_dict, 'workflow/scripts/wt.fasta', int(epi_lo), int(epi_up))
    # mt_binding_affinity = calc_peptide_binding(alleles_dict, 'workflow/scripts/mt.fasta', int(epi_lo), int(epi_up))

    # for x in out_itrm:
        # if x['wt_pep'] in wt_binding_affinity:
            # print(x)
            # print(wt_binding_affinity[x['wt_pep']])


    






#                                for x in alleles_dict.keys():
                                    # group = x
                                    # alleles_per_group = alleles_dict[x]
                                    # for allele in alleles_per_group:
                                        # allele_mod = 'HLA-'+allele
                                        # if valid_alleles(allele_mod):
                                            # if '$' not in subseqs[key][0]:
                                                # wt_ic50, wt_rank = calc_peptide_binding(allele_mod, subseqs[key][0])
                                                # wt_immuno = calc_immunogenecity(subseqs[key][0])
                                            # else:
                                                # wt_ic50 = wt_rank = wt_immuno = 'NA'

                                            # mt_ic50, mt_rank = calc_peptide_binding(allele_mod, subseqs[key][1])
                                            # mt_immuno = calc_immunogenecity(subseqs[key][1])

                                            # print(mt_ic50)
                                            # print(mt_immuno)
                                            # if float(mt_ic50) <= 500:
                                                # output_row['wt_sequence'] = subseqs[key][0]
                                                # output_row['mt_sequence'] = subseqs[key][1]
                                                # output_row['wt_ic50'] = wt_ic50
                                                # output_row['mt_ic50'] = mt_ic50
                                                # output_row['wt_immuno'] = allele_mod
                                                # output_row['mt_immuno'] = allele_mod
                                                # output_row['allele'] = allele_mod

                                                # print_entries(output_row, outputfile)
