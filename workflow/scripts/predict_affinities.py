import sys
import os
import subprocess
import configargparse
import tempfile

"""
Usage: python predict_affinities.py -i <peptides> \
        -a <alleles> \
        -e <epitope_lengths> \
        -o <output> \
        -l <log>
"""

def main():
    options = parse_arguments()

    with tempfile.TemporaryDirectory() as tmp_seqs:
        alleles = []
        fh_alleles = open(options.alleles, 'r')
        for allele in fh_alleles:
            alleles.append(allele.rstrip())

        wt_fname = {}
        mt_fname = {}

        # counters to store number of sequence in fa
        wt_cnt = {}
        mt_cnt = {}
        # file handler for wt and mt sequences
        fh_wt = {}
        fh_mt = {}

        epilens = options.lengths.split(',')
        epilens = [eval(i) for i in epilens]
        for epilen in epilens:

            wt_fname[epilen] = os.path.join(tmp_seqs, f'wt_{epilen}.fa')
            mt_fname[epilen] = os.path.join(tmp_seqs, f'mt_{epilen}.fa')

            fh_wt[epilen] = open(wt_fname[epilen], 'w')
            fh_mt[epilen] = open(mt_fname[epilen], 'w')

            # initialise counter
            wt_cnt[epilen] = 1
            mt_cnt[epilen] = 1

        subseqs = []

        # iterate through file
        fh = open(options.input, 'r')
        next(fh)   # skip header
        for line in fh:
            entries = line.rstrip().split('\t')
            # make sure that 
            if '$' in entries[11]:
                wt_epitope_seq = entries[11].split('$')[0]
            else:
                wt_epitope_seq = entries[11]
            mt_epitope_seq = entries[12]

            for epilen in epilens:
                if len(wt_epitope_seq) >= epilen+1:
                    fh_wt[epilen].write(f'>{wt_cnt[epilen]}\n{wt_epitope_seq}\n')
                    entries.append(wt_cnt[epilen])
                    wt_cnt[epilen] += 1
                else:
                    entries.append(0)

                if len(mt_epitope_seq) >= epilen+1:
                    fh_mt[epilen].write(f'>{mt_cnt[epilen]}\n{mt_epitope_seq}\n')
                    entries.append(mt_cnt[epilen])
                    mt_cnt[epilen] += 1
                else:
                    entries.append(0)

            subseqs.append(entries)

        # its necessary to close the files
        fh.close()
        for epilen in epilens:
            fh_wt[epilen].close()
            fh_mt[epilen].close()
        
        wt_affinities = calc_peptide_binding(alleles, wt_fname, 'wt')
        mt_affinities = calc_peptide_binding(alleles, mt_fname, 'mt')

#        print(mt_affinities)

#        print(wt_affinities.keys())

        outputfile = open(options.output, 'w')
        print_header(outputfile)

        for entry in subseqs:
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
            final_result['vaf'] = float(entry[14])
            final_result['wt_allele_depth'] = int(entry[15])
            final_result['mt_allele_depth'] = int(entry[16])
            final_result['read_depth'] = int(entry[17])

            # check results of different epilens
            for epilen_idx in range(0,len(epilens)):

                # retrieve sequence numbers (from binding affinities results)
                wt_seqnum = int(entry[18:][epilen_idx*2])
                mt_seqnum = int(entry[18:][epilen_idx*2+1])

                wt = None
                if wt_seqnum in wt_affinities[epilens[epilen_idx]].keys():
                    wt = wt_affinities[epilens[epilen_idx]][wt_seqnum]
                else:
                    final_result['wt_epitope_ic50'] = 'NA'
                    final_result['wt_epitope_rank'] = 'NA'

                if mt_seqnum in mt_affinities[epilens[epilen_idx]].keys():
                    mt = mt_affinities[epilens[epilen_idx]][mt_seqnum]
                else:
                    continue

                for epitope in mt.keys(): 
                    # determine by mhc_i / convert to 0-based
                    start_pos_in_subseq = int(mt[epitope][1])
                    end_pos_in_subseq = int(mt[epitope][2])

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

                    # # search for corresponding WT
                    startpos_epitope_subseq = final_result['mt_epitope_seq'].find(epitope)
                    final_result['wt_epitope_seq'] = final_result['wt_subseq'][startpos_epitope_subseq:len(epitope)]

                    if wt is not None:
                        if final_result['wt_epitope_seq'] in wt.keys():
                            final_result['wt_epitope_ic50'] = wt[final_result['wt_epitope_seq']][3]
                            final_result['wt_epitope_rank'] = wt[final_result['wt_epitope_seq']][4]
                        else:
                            final_result['wt_epitope_ic50'] = 'NA'
                            final_result['wt_epitope_rank'] = 'NA'


                    # calculate ranking calc_ranking_score 
                    score = calc_ranking_score(final_result['vaf'], final_result['wt_epitope_ic50'], final_result['mt_epitope_ic50'])
                    final_result['ranking_score'] = score

                    print_row(final_result, outputfile)



def calc_peptide_binding(alleles, fnames, group):
    binding_affinities = {}
    for allele in alleles:
        for epilen in fnames.keys():
            print(f'calculate binding affinities - allele:{allele} epilen:{epilen} group: {group}')
            binding_affinities[epilen] = {}
            call = ['python', 
                'workflow/scripts/mhc_i/src/predict_binding.py', 
                'netmhcpan', 
                allele, 
                str(epilen),
                fnames[epilen]]

            result = subprocess.run(call,
                stdout = subprocess.PIPE,
                universal_newlines=True)
            predictions = result.stdout.rstrip().split('\n')[1:]
            for line in predictions:
                entries = line.split('\t')
                if group == 'mt':
                    if float(entries[8]) >= 500:
                        continue

                # sequence number in results used as key
                seqnum = int(entries[1])
                epitope_seq = entries[5]
                allele = entries[0]

                # start and end in sequence (0-based)
                start = int(entries[2])-1 
                end = int(entries[3])-1

                ic50 = float(entries[8])
                rank = float(entries[9])

                if seqnum not in binding_affinities[epilen]:
                    binding_affinities[epilen][seqnum] = {}

                if entries[5] not in binding_affinities[epilen][seqnum]:
                    binding_affinities[epilen][seqnum][epitope_seq] = (allele, start, end, ic50, rank)

    return binding_affinities



def calc_ranking_score(vaf, wt_ic50, mt_ic50):
    # fold change
    if wt_ic50 != 'NA':
        fold_change = float(wt_ic50) / float(mt_ic50)
    else:
        fold_change = float(100000.0) / float(mt_ic50)

    score = (1/float(mt_ic50)) + fold_change + vaf * 100.0

    return score

            
def print_header(outputfile):
    outputfile.write(f"chrom\tstart\tend\tref\talt\tgene_name\tgene_id")
    outputfile.write(f"\ttranscript_id\tsource\tgroup\tvariant_type")
    outputfile.write(f"\tallele\twt_epitope_seq\twt_epitope_ic50\twt_epitope_rank")
    outputfile.write(f"\tmt_epitope_seq\tmt_epitope_ic50\tmt_epitope_rank\tmutation_position")
    outputfile.write(f"\tvaf\twt_allele_depth\tmt_allele_depth\tread_depth\t")
    outputfile.write(f"\tranking_score\n")


def print_row(row, outputfile):
    outputfile.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['ref']}\t{row['alt']}\t{row['gene_name']}\t{row['gene_id']}")
    outputfile.write(f"\t{row['transcript_id']}\t{row['source']}\t{row['group']}\t{row['variant_type']}")
    outputfile.write(f"\t{row['allele']}\t{row['wt_epitope_seq']}\t{row['wt_epitope_ic50']}\t{row['wt_epitope_rank']}")
    outputfile.write(f"\t{row['mt_epitope_seq']}\t{row['mt_epitope_ic50']}\t{row['mt_epitope_rank']}\t{row['mutation_position']}")
    outputfile.write(f"\t{row['vaf']}\t{row['wt_allele_depth']}\t{row['mt_allele_depth']}\t{row['read_depth']}")
    outputfile.write(f"\t{row['ranking_score']}\n")

def parse_arguments():
    p = configargparse.ArgParser()
    p.add('-i', '--input', required=True, help='Input vcf files')
    p.add('-a', '--alleles', required=True, help='Alleles to be predicted')
    p.add('-e', '--lengths', required=True, help='Lengths of epitopes to be predicted')
    p.add('-o', '--output', required=True, help='Output table')
    p.add('-l', '--log', required=False, help='Log file')

    options = p.parse_args()

    return options

main()
