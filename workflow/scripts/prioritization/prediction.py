import tempfile
import os
import concurrent.futures
import subprocess
from pathlib import Path

class BindingAffinities:
    def __init__(self, threads):
        self.threads = threads
        pass

    def start(self, allele_file, epitope_lengths, output_dir, mhc_class):
        # create temorary_directory
        with tempfile.TemporaryDirectory() as tmp_seqs:
            self.get_alleles(allele_file)

            # initialize filenames
            wt_fname = {}
            mt_fname = {}

            # counters to store number of sequence in fa
            wt_cnt = {}
            mt_cnt = {}

            # file handler for wt and mt sequences
            fh_wt = {}
            fh_mt = {}

            epilens = self.extract_epilens(epitope_lengths)

            for epilen in epilens:
                wt_fname[epilen] = os.path.join(tmp_seqs, f'wt_{epilen}.fa')
                mt_fname[epilen] = os.path.join(tmp_seqs, f'mt_{epilen}.fa')

                fh_wt[epilen] = open(wt_fname[epilen], 'w')
                fh_mt[epilen] = open(mt_fname[epilen], 'w')

                # initialise counter
                wt_cnt[epilen] = 1
                mt_cnt[epilen] = 1
            
            subseqs = []

            # iterate through the varianteffectsfile
            fh = open(Path(output_dir, "variant_effects.tsv"), 'r')
            next(fh)   # skip header
            for line in fh:
                entries = line.rstrip().split('\t')
                for epilen in epilens:
                    aa_var_start = int(entries[12])
                    aa_var_end = int(entries[13])

                    wt_subseq = entries[9]
                    mt_subseq = entries[10]

                    # adjust the length of the subsequence according to epilen
                    if aa_var_start >= epilen + 1:
                        left = aa_var_start - (epilen - 1)
                    else:
                        left = 0

                    if aa_var_end + (epilen - 1) <= len(mt_subseq):
                        right = aa_var_end + (epilen - 1)
                    else:
                        right = len(mt_subseq) - 1 

                    """ IEDB requires the sequence to be at least the length
                    of the epitope + 1 to function - if this is not satisfied
                    try to extend it to the left """
                    while left > 0 and right - left + 1 < epilen + 1:
                        left -= 1
#                    if right - left + 1 < epilen + 1:
#                        continue

                    wt_subseq_adj = wt_subseq[left:right+1]
                    mt_subseq_adj = mt_subseq[left:right+1]

                    # determine the epitope sequences
                    if '$' in wt_subseq_adj:
                        wt_epitope_seq = wt_subseq_adj.split('$')[0]
                    else:
                        wt_epitope_seq = wt_subseq_adj
                    mt_epitope_seq = mt_subseq_adj

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

                    #print(entries)

                subseqs.append(entries)

            # its necessary to close the files
            fh.close() # variant_effects
            for epilen in epilens:
                fh_wt[epilen].close()
                fh_mt[epilen].close()

            print(f"calculate binding affinities...")
            wt_affinities = self.collect_binding_affinities(self.alleles, 
                                                            wt_fname, 
                                                            epilens, 
                                                            'wt',
                                                            mhc_class,
                                                            self.threads)
            
            mt_affinities = self.collect_binding_affinities(self.alleles,
                                                            mt_fname,
                                                            epilens,
                                                            'mt',
                                                            mhc_class,
                                                            self.threads)
            print("Done")
            
            outfile = open(os.path.join(output_dir,
                                        f'{mhc_class}_neoepitopes.txt'),'w')
            BindingAffinities.write_header(outfile)


            for entry in subseqs:
                final = {}
                final["chrom"] = entry[0]
                final["start"] = entry[1]
                final["end"] = entry[2]
                final["gene_id"] = entry[3]
                final["gene_name"] = entry[4]
                final["transcript_id"] = entry[5]
                final["source"] = entry[6]
                final["group"] = entry[7]
                final["var_type"] = entry[8]
                final["var_start"] = entry[11]
                final["vaf"] = float(entry[14])
                final["supp"] = entry[15]
                final["depth"] = entry[16]
                final["TPM"] = entry[17]
                final["NMD"] = entry[18]
                final["PTC_dist_ejc"] = entry[19]
                final["PTC_exon_number"] = entry[20]
                final["NMD_escape_rule"] = entry[21]

                aa_var_start = int(entry[12])
                aa_var_end = int(entry[13])

                # extract the subsequences (needed to determine the wt epitope)
                wt_subseq = entry[9]
                mt_subseq = entry[10]

                for epilen_idx in range(0, len(epilens)):

                    # the sequence number of each entry corresponds to the 
                    # sequence number of the epitope in the fasta file
                    wt_seqnum = int(entry[22:][epilen_idx*2])
                    mt_seqnum = int(entry[22:][epilen_idx*2+1])

                    print(f"epilen_idx:{epilen_idx}")
                    print(wt_affinities)




                    wt = None
                    if wt_seqnum in wt_affinities[epilens[epilen_idx]].keys():
                        wt = wt_affinities[epilens[epilen_idx]][wt_seqnum]
                    else:
                        final["wt_epitope_ic50"] = None
                        final["wt_epitope_rank"] = None

                    if mt_seqnum in mt_affinities[epilens[epilen_idx]].keys():
                        mt = mt_affinities[epilens[epilen_idx]][mt_seqnum]
                    else:
                        continue

                    for epitope in mt.keys():
                        # (allele, start, end, ic50, rank)
                        # determine by mhc_i (as determined by IEDB - 0 based)
                        start_pos_in_subseq = int(mt[epitope][1])
                        end_pos_in_subseq = int(mt[epitope][2])

                        # check if the mutation (aa_var_[start|end]) is either 
                        # part of the epitope or upstream of it - this ensures 
                        # that the sequence is mutated
                        # if (aa_var_end < start_pos_in_subseq or
                            # aa_var_start > end_pos_in_subseq):
                            # continue

                        final["mt_epitope_seq"] = epitope
                        final["allele"] = mt[epitope][0]
                        final["mt_epitope_ic50"] = mt[epitope][3]
                        final["mt_epitope_rank"] = mt[epitope][4]

                        """ search for corresponding WT by first searching for 
                        the epitope in the mt_subseq and using this information
                        to return the WT epitope sequence"""
                        startpos_epitope_in_subseq = mt_subseq.find(epitope)
                        startpos = startpos_epitope_in_subseq
                        final["wt_epitope_seq"] = wt_subseq[startpos:startpos+len(epitope)]

                        # search for binidng affinities of wildtype sequence
                        final["wt_epitope_ic50"] = None
                        final["wt_epitope_rank"] = None
                        if wt is not None:
                            if final["wt_epitope_seq"] in wt.keys():
                                final["wt_epitope_ic50"] = wt[final["wt_epitope_seq"]][3]
                                final["wt_epitope_rank"] = wt[final["wt_epitope_seq"]][4]

                        # calculate ranking calc_ranking_score 
                        score = BindingAffinities.calc_ranking_score(final['vaf'], 
                                                                     final['wt_epitope_ic50'], 
                                                                     final['mt_epitope_ic50'])
                        final['ranking_score'] = score
                        final['agretopicity'] = BindingAffinities.calc_agretopicity(final["wt_epitope_ic50"],
                                                                                    final["mt_epitope_ic50"])

                    
                        BindingAffinities.write_entry(final, outfile)
            outfile.close()



    def get_alleles(self, allele_file):
        self.alleles = {}
        fh_alleles = open(allele_file, 'r')
        for allele in fh_alleles:
            line = allele.rstrip().split('\t')
            self.alleles[line[1]] = line[0]

    @staticmethod
    def extract_epilens(lengths):
        """determines the lengths of the epitopes specified of the commandline
        as these can be specific a single digits (e.g, 8,9,10) or ranges 
        (e.g., 8-10) this function extracts the single lengths"""
        epilens = []
        for lens in lengths.split(','):
            if '-' in lens:
                lower, upper = lens.split('-', 1)
                epilens.extend(range(int(lower), int(upper)+1))
            else:
                epilens.append(int(lens))
        return epilens
    
    @staticmethod
    def collect_binding_affinities(alleles, 
                                   fnames, 
                                   epilens, 
                                   group, 
                                   mhc_class, 
                                   threads):
        affinities_results = {}
        number_threads = int(threads)

        # # only iterate through keys (alleles)
        # # TODO: also incorporate sources
        with concurrent.futures.ThreadPoolExecutor(max_workers=number_threads) as executor:
            futures = {}
            for allele in alleles:
                for epilen in epilens:
                    future = executor.submit(BindingAffinities.calc_binding_affinities, 
                                             fnames[epilen], 
                                             allele, 
                                             epilen, 
                                             group,
                                             mhc_class)
                    futures[future] = epilen


        for future in concurrent.futures.as_completed(futures):
            epilen = futures[future]
            binding_affinities = future.result()

            # check epilen already present
            if epilen not in affinities_results.keys():
                affinities_results[epilen] = {}

            for seqnum in binding_affinities.keys():
                if seqnum not in affinities_results[epilen].keys():
                    affinities_results[epilen][seqnum] = binding_affinities[seqnum]
                else:
                    for seq in binding_affinities[seqnum].keys():
                        if seq not in affinities_results[epilen][seqnum].keys():
                            affinities_results[epilen][seqnum][seq] = binding_affinities[seqnum][seq]

        return affinities_results
    
    @staticmethod
    def calc_binding_affinities(fa_file, allele, epilen, group, mhc_class):
        binding_affinities = {}
# #    binding_affinities[epilen] = {}
        print(f'calculate binding affinities - allele:{allele} epilen:{epilen} group: {group}')
        call = ['python']
        if mhc_class == "mhc-I":
            call.append("workflow/scripts/mhc_i/src/predict_binding.py")
            call.append("netmhcpan")
        elif mhc_class == "mhc-II":
            call.append("workflow/scripts/mhc_ii/mhc_II_binding.py")
            call.append("netmhciipan_ba")
        call.append(allele)
        call.append(str(epilen))
        call.append(fa_file)

        # make sure that there are entries in the file
        if os.stat(fa_file).st_size != 0:
            result = subprocess.run(call,
                stdout = subprocess.PIPE,
                universal_newlines=True)
            predictions = result.stdout.rstrip().split('\n')[1:]

            for line in predictions:
                entries = line.split('\t')
                if group == 'mt':
                    if float(entries[8]) >= 500:
                        continue

                # start and end in sequence (0-based)
                start = int(entries[2])-1 
                end = int(entries[3])-1
                
                # sequence number in results used as key
                allele = entries[0]
                seqnum = int(entries[1])
                if mhc_class == "mhc-I":
                    epitope_seq = entries[5]
                    ic50 = float(entries[8])
                    rank = float(entries[9])
                elif mhc_class == "mhc-II":
                    epitope_seq = entries[6]
                    ic50 = float(entries[7])
                    rank = float(entries[8])

                if seqnum not in binding_affinities:
                    binding_affinities[seqnum] = {}

                if epitope_seq not in binding_affinities[seqnum]:
                    binding_affinities[seqnum][epitope_seq] = (allele, start, end, ic50, rank)

        return binding_affinities
    
    @staticmethod
    def calc_ranking_score(vaf, wt_ic50, mt_ic50):
        if vaf == -1.0:
            return -1.0

        # fold change
        if wt_ic50 != None:
            fold_change = float(mt_ic50) / float(wt_ic50)
        else:
            fold_change = float(100000.0) / float(mt_ic50)

        score = (1/float(mt_ic50)) + fold_change + vaf * 100.0

        return score


    def calc_agretopicity(wt_ic50, mt_ic50):
        if wt_ic50 == None:
            return None
        else:
            return float(mt_ic50) / float(wt_ic50)


    @staticmethod
    def write_header(outputfile):
        outputfile.write(f"chrom\tstart\tend\tallele\tgene_id\tgene_name\t")
        outputfile.write(f"transcript_id\tsource\tgroup\tvar_type\t")
        outputfile.write(f"var_start\twt_epitope_seq\twt_epitope_ic50\t")
        outputfile.write(f"wt_epitope_rank\tmt_epitope_seq\t")
        outputfile.write(f"mt_epitope_ic50\tmt_epitope_rank\tvaf\tsupporting\t")
        outputfile.write(f"TPM\tagretopicity\tNMD\tPTC_dist_ejc\t")
        outputfile.write(f"PTC_exon_number\tNMD_escape_rule\n")

    @staticmethod
    def write_entry(row, outputfile):
        outputfile.write(f'{format_output(row["chrom"])}\t')
        outputfile.write(f'{format_output(row["start"])}\t')
        outputfile.write(f'{format_output(row["end"])}\t')
        outputfile.write(f'{format_output(row["allele"])}\t')
        outputfile.write(f'{format_output(row["gene_id"])}\t')
        outputfile.write(f'{format_output(row["gene_name"])}\t')
        outputfile.write(f'{format_output(row["transcript_id"])}\t')
        outputfile.write(f'{format_output(row["source"])}\t')
        outputfile.write(f'{format_output(row["group"])}\t')
        outputfile.write(f'{format_output(row["var_type"])}\t')
        outputfile.write(f'{format_output(row["var_start"])}\t')
        outputfile.write(f'{format_output(row["wt_epitope_seq"])}\t')
        outputfile.write(f'{format_output(row["wt_epitope_ic50"])}\t')
        outputfile.write(f'{format_output(row["wt_epitope_rank"])}\t')
        outputfile.write(f'{format_output(row["mt_epitope_seq"])}\t')
        outputfile.write(f'{format_output(row["mt_epitope_ic50"])}\t')
        outputfile.write(f'{format_output(row["mt_epitope_rank"])}\t')
        outputfile.write(f'{format_output(row["vaf"])}\t')
        outputfile.write(f'{format_output(row["supp"])}\t')
        outputfile.write(f'{format_output(row["TPM"])}\t')
        outputfile.write(f'{format_output(row["agretopicity"])}\t')
        outputfile.write(f'{format_output(row["NMD"])}\t')
        outputfile.write(f'{format_output(row["PTC_dist_ejc"])}\t')
        outputfile.write(f'{format_output(row["PTC_exon_number"])}\t')
        outputfile.write(f'{format_output(row["NMD_escape_rule"])}\n')
        
    
def format_output(field):
    if field == None:
        return '.'
    else:
        return field
