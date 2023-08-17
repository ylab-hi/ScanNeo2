import sys
import subprocess

"""
Usage: python predict_affinities.py <peptides> <alleles> <epitope_len>
"""


def main():
    
    peptides = []
    #generate fasta
    wt_fasta = open('workflow/scripts/wt.fa','w')
    mt_fasta = open('workflow/scripts/mt.fa','w')

    wt_seqnum = 1
    mt_seqnum = 1

    # laod peptides
    fh = open(sys.argv[1], 'r')
    next(fh)
    for line in fh:
        entries = line.rstrip().split('\t')
        
        if '$' in entries[11]:
            wt_seq = entries[11].split('$')[0]
        else:
            wt_seq = entries[11]

        if len(wt_seq) >= 9:
            wt_fasta.write(f'>{wt_seqnum}\n{wt_seq}\n')
            entries.append(wt_seqnum)
            wt_seqnum += 1
        else:
            entries.append(0)

        mt_seq = entries[12]
        if len(mt_seq) >= 9:
            mt_fasta.write(f'>{mt_seqnum}\n{mt_seq}\n')
            entries.append(mt_seqnum)
            mt_seqnum += 1
        else:
            entries.append(0)


        peptides.append(entries)

    fh.close()
    wt_fasta.close()
    mt_fasta.close()

    # parse alleles from input
    alleles = []
    al_handler = open(sys.argv[2], 'r')
    for al in al_handler:
        alleles.append(al.rstrip())

    # parse epitope length
    epilens = sys.argv[3].split(',')

    wt_affinities = calc_peptide_binding(alleles, 'workflow/scripts/wt.fa', epilens, 'wt')
    mt_affinities = calc_peptide_binding(alleles, 'workflow/scripts/mt.fa', epilens, 'mt')

    outputfile = open(sys.argv[4], 'w')
    print_header(outputfile)

    for entry in peptides:
        final_result = {}
        final_result['chrom'] = entry[0]
        final_result['start'] = int(entry[1])
        final_result['end'] = int(entry[2])
        final_result['ref'] = entry[3]
        final_result['alt'] = entry[4]
        final_result['gene_name'] = entry[5]
        final_result['gene_id'] = entry[6]
        final_result['transcript_id'] = entry[7]
        final_result['source'] = entry[8]
        final_result['group'] = entry[9]
        final_result['variant_type'] = entry[10]
        final_result['wt_subseq'] = entry[11]
        final_result['mt_subseq'] = entry[12]
        final_result['new_var_start_pos'] = int(entry[13])
        final_result['wt_seqnum'] = entry[14]
        final_result['mt_seqnum'] = entry[15]

        wt = None
        # search for epitopes with high binding affinity
        if final_result['wt_seqnum'] in wt_affinities:
            wt = wt_affinities[final_result['wt_seqnum']]
        else:
            final_result['wt_epitope_ic50'] = 'NA'
            final_result['wt_epitope_rank'] = 'NA'
        
        if final_result['mt_seqnum'] in mt_affinities:
            mt = mt_affinities[final_result['mt_seqnum']]
        else:
            continue

        for epitope in mt.keys():
            # determine by mhc_i / convert to 0-based
            start_pos_in_subseq = int(mt[epitope][1]-1)
            end_pos_in_subseq = int(mt[epitope][2]-1)

            # check if mutation is part of the subsequence (within or upstream)
            if final_result['new_var_start_pos'] >= end_pos_in_subseq:
                continue
            elif final_result['new_var_start_pos'] <= start_pos_in_subseq:
                final_result['mutation_position'] = 0
            else:
                final_result['mutation_position'] = final_result['new_var_start_pos'] - start_pos_in_subseq

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

            print_row(final_result, outputfile)


def print_header(outputfile):
    outputfile.write(f"chrom\tstart\tend\tref\talt\tgene_name\tgene_id")
    outputfile.write(f"\ttranscript_id\tsource\tgroup\tvariant_type")
    outputfile.write(f"\tallele\twt_epitope_seq\twt_epitope_ic50\twt_epitope_rank")
    outputfile.write(f"\tmt_epitope_seq\tmt_epitope_ic50\tmt_epitope_rank\tmutation_position\n")

def print_row(row, outputfile):
    outputfile.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['ref']}\t{row['alt']}\t{row['gene_name']}\t{row['gene_id']}")
    outputfile.write(f"\t{row['transcript_id']}\t{row['source']}\t{row['group']}\t{row['variant_type']}")
    outputfile.write(f"\t{row['allele']}\t{row['wt_epitope_seq']}\t{row['wt_epitope_ic50']}\t{row['wt_epitope_rank']}")
    outputfile.write(f"\t{row['mt_epitope_seq']}\t{row['mt_epitope_ic50']}\t{row['mt_epitope_rank']}\t{row['mutation_position']}\n")

def calc_peptide_binding(alleles, fa_file, epilens, wt_mt):
    binding_affinity = {}
    for allele in alleles:
        for epilen in epilens:
            print(f"allele: {allele}, epilen: {epilen}")
            call = ['python', 
                'workflow/scripts/mhc_i/src/predict_binding.py', 
                'netmhcpan', 
                allele, 
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

main()
