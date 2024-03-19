import tempfile
import os
import concurrent.futures
import subprocess
from pathlib import Path
import utility as ut


class PredictionResults:
    def __init__(self, epilens, tmp):
        self.tmp = tmp
        self.max_entries = 100
        self.fnames, self.counter = self.initialize(epilens)

    def add_pepseq(self, pepseq, epilen):
        if self.counter[epilen][1] < self.max_entries:
            self.add_fasta_entry(pepseq, epilen)
            self.counter[epilen][1] += 1 # increase entry counter

        elif self.counter[epilen][1] == self.max_entries: # file is full
            self.counter[epilen][0] += 1 # go to the next file
            self.counter[epilen][1] = 1 # reset entry counter
            tmpfile = os.path.join(self.tmp, f"{self.counter[epilen][0]}_{epilen}.fa")
            self.fnames[epilen][self.counter[epilen][0]] = tmpfile
            self.add_fasta_entry(pepseq, epilen)

        return self.counter[epilen][0], self.counter[epilen][1]

    def add_fasta_entry(self, pepseq, epilen):
        fh = open(self.fnames[epilen][self.counter[epilen][0]], 'a')
        fh.write(f">{self.counter[epilen][1]}\n{pepseq}\n")
        fh.close()

    # initalize class variables/files
    def initialize(self, epilens):
        fnames = {}
        counter = {}
        for epilen in epilens:
            # counter (file and entry)
            counter[epilen] = [0,0] # file counter, entry 
            fnames[epilen] = {}
            tmpfile = os.path.join(self.tmp, f"{counter[epilen][0]}_{epilen}.fa")
            fnames[epilen][counter[epilen][0]] = tmpfile

        return fnames, counter


class BindingAffinities:
    def __init__(self, threads):
        self.threads = threads
        pass

    def start(self, allele_file, epitope_lengths, output_dir, mhc_class, vartype):
        # create temorary_directory
        with tempfile.TemporaryDirectory() as tmp_seqs:
            alleles = ut.get_alleles(allele_file)

            # number of entries in the fasta file (to distribute for prediction)
            number_entries = 100

            epilens = ut.extract_epilens(epitope_lengths)
            wt_pred = PredictionResults(epilens, tmp_seqs)
            mt_pred = PredictionResults(epilens, tmp_seqs)

            subseqs = []

            # iterate through the varianteffectsfile
            fh = open(Path(output_dir, f"{vartype}_variant_effects.tsv"), 'r')
            next(fh)   # skip header
            for line in fh:
                entries = line.rstrip().split('\t')
                    
                aa_var_start = int(entries[12])
                aa_var_end = int(entries[13])

                wt_subseq = entries[9]
                mt_subseq = entries[10]

                for epilen in epilens:
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
                    try to extend it to the left (while still including the mutation) """
                    while left > 0 and right - left + 1 < epilen + 1:
                        left -= 1
                    
                    wt_subseq_adj = wt_subseq[left:right+1]
                    mt_subseq_adj = mt_subseq[left:right+1]
                    
                    # determine the epitope sequences
                    if '$' in wt_subseq_adj:
                        wt_epitope_seq = wt_subseq_adj.split('$')[0]
                    else:
                        wt_epitope_seq = wt_subseq_adj
                    mt_epitope_seq = mt_subseq_adj
                    
#                    print(f"wt_epitope_seq: {wt_epitope_seq}")
#                    print(f"mt_epitope_seq: {mt_epitope_seq}")
                    
                    if len(wt_epitope_seq) >= epilen+1:
                        fcnt, ecnt = wt_pred.add_pepseq(wt_epitope_seq, epilen)
                        entries.append(f"{fcnt}_{ecnt}")
                    else:
                        entries.append(f"{-1}_{-1}")
                    
                    if len(mt_epitope_seq) >= epilen+1:
                        fcnt, ecnt = mt_pred.add_pepseq(mt_epitope_seq, epilen)
                        entries.append(f"{fcnt}_{ecnt}")
                    else:
                        entries.append(f"{-1}_{-1}")

                subseqs.append(entries)
            
            print(f"calculate binding affinities...")
            wt_affinities = self.collect_binding_affinities(alleles, 
                                                            wt_pred.fnames,
                                                            epilens,
                                                            'wt',
                                                            mhc_class,
                                                            self.threads)
            
            mt_affinities = self.collect_binding_affinities(alleles, 
                                                            mt_pred.fnames,
                                                            epilens,
                                                            'mt',
                                                            mhc_class,
                                                            self.threads)
            print("...done")
            
            outfile = open(os.path.join(output_dir,
                                        f"{vartype}_{mhc_class}_neoepitopes.txt"),"w")
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
                    wt_nkey = int(entry[22:][epilen_idx*2])
                    mt_nkey = int(entry[22:][epilen_idx*2+1])

                    wt = None
                    if wt_nkey in wt_affinities[epilens[epilen_idx]].keys():
                        wt = wt_affinities[epilens[epilen_idx]][wt_nkey]
                    else:
                        final["wt_epitope_ic50"] = None
                        final["wt_epitope_rank"] = None

                    if mt_nkey in mt_affinities[epilens[epilen_idx]].keys():
                        mt = mt_affinities[epilens[epilen_idx]][mt_nkey]
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
                    for fkey in fnames[epilen].keys():
                        future = executor.submit(BindingAffinities.calc_binding_affinities,
                                                 fnames[epilen][fkey],
                                                 fkey,
                                                 allele,
                                                 epilen,
                                                 group,
                                                 mhc_class)
                        futures[future] = [fkey, epilen]
        
        for future in concurrent.futures.as_completed(futures):
            fkey, epilen = futures[future]
            binding_affinities = future.result()

            # check if epilen already present
            if epilen not in affinities_results.keys():
                affinities_results[epilen] = {}

            for seqnum in binding_affinities.keys():
                nkey = f"{fkey}_{seqnum}"
                if nkey not in affinities_results[epilen].keys():
                    affinities_results[epilen][nkey] = binding_affinities[seqnum]
                else:
                    for seq in binding_affinities[seqnum].keys():
                        if seq not in affinities_results[epilen][nkey].keys():
                            affinities_results[epilen][nkey][seq] = binding_affinities[seqnum][seq]

        return affinities_results

    
    @staticmethod
    def calc_binding_affinities(fa_file, fkey, allele, epilen, group, mhc_class):
        binding_affinities = {}
        print(f"Calculate binding affinities - allele: {allele} epilen: {epilen} group: {group} file: {fkey}") 
        call = ["python"]
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
        outputfile.write(f'{ut.format_output(row["chrom"])}\t')
        outputfile.write(f'{ut.format_output(row["start"])}\t')
        outputfile.write(f'{ut.format_output(row["end"])}\t')
        outputfile.write(f'{ut.format_output(row["allele"])}\t')
        outputfile.write(f'{ut.format_output(row["gene_id"])}\t')
        outputfile.write(f'{ut.format_output(row["gene_name"])}\t')
        outputfile.write(f'{ut.format_output(row["transcript_id"])}\t')
        outputfile.write(f'{ut.format_output(row["source"])}\t')
        outputfile.write(f'{ut.format_output(row["group"])}\t')
        outputfile.write(f'{ut.format_output(row["var_type"])}\t')
        outputfile.write(f'{ut.format_output(row["var_start"])}\t')
        outputfile.write(f'{ut.format_output(row["wt_epitope_seq"])}\t')
        outputfile.write(f'{ut.format_output(row["wt_epitope_ic50"])}\t')
        outputfile.write(f'{ut.format_output(row["wt_epitope_rank"])}\t')
        outputfile.write(f'{ut.format_output(row["mt_epitope_seq"])}\t')
        outputfile.write(f'{ut.format_output(row["mt_epitope_ic50"])}\t')
        outputfile.write(f'{ut.format_output(row["mt_epitope_rank"])}\t')
        outputfile.write(f'{ut.format_output(row["vaf"])}\t')
        outputfile.write(f'{ut.format_output(row["supp"])}\t')
        outputfile.write(f'{ut.format_output(row["TPM"])}\t')
        outputfile.write(f'{ut.format_output(row["agretopicity"])}\t')
        outputfile.write(f'{ut.format_output(row["NMD"])}\t')
        outputfile.write(f'{ut.format_output(row["PTC_dist_ejc"])}\t')
        outputfile.write(f'{ut.format_output(row["PTC_exon_number"])}\t')
        outputfile.write(f'{ut.format_output(row["NMD_escape_rule"])}\n')
        
    
