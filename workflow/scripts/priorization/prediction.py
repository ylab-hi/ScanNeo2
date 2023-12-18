import tempfile
import pdb
import os
import concurrent.futures
import subprocess

class BindingAffinities:
    def __init__(self, effects, mhcI, mhcII, mhcI_len, mhcII_len, threads, outdir):
    
        with tempfile.TemporaryDirectory() as tmp_seqs:
            if mhcI is not None:
                mhcI_alleles = {}
                fh_mhcI_alleles = open(mhcI, 'r')
                for allele in fh_mhcI_alleles:
                    line = allele.rstrip().split('\t')
                    mhcI_alleles[line[1]] = line[0]

            wt_fname = {}
            mt_fname = {}
        
            # counters to store number of sequence in fa
            wt_cnt = {}
            mt_cnt = {}

            # file handler for wt and mt sequences
            fh_wt = {}
            fh_mt = {}

            # extract mhcI lens
            epilens = []
            for lens in mhcI_len.split(','):
                if '-' in lens:
                    lower, upper = lens.split('-',1)
                    epilens.extend(range(int(lower), int(upper)+1))
                else:
                    epilens.append(int(lens))

            for epilen in epilens:
                wt_fname[epilen] = os.path.join(tmp_seqs, f'wt_{epilen}.fa')
                mt_fname[epilen] = os.path.join(tmp_seqs, f'mt_{epilen}.fa')

                fh_wt[epilen] = open(wt_fname[epilen], 'w')
                fh_mt[epilen] = open(mt_fname[epilen], 'w')

                # initialise counter
                wt_cnt[epilen] = 1
                mt_cnt[epilen] = 1

            subseqs = []

            # iterate through the file
            fh = open(effects, 'r')
            next(fh)   # skip header
            for line in fh:
                entries = line.rstrip().split('\t')
                # print(line)


                for epilen in epilens:
                    local_var_start = int(entries[15])
                    local_var_end = int(entries[16])

                    wt_subseq = entries[12]
                    mt_subseq = entries[13]

                    # adjust length of seqs (according to epilen)
                    if local_var_start >= epilen + 1:
                        left = local_var_start - (epilen - 1)
                    else:
                        left = 0

                    if local_var_end + (epilen - 1) <= len(mt_subseq):
                        right = local_var_end + epilen - 1
                    else:
                        right = len(mt_subseq)

                    wt_subseq_adj = wt_subseq[left:right]
                    mt_subseq_adj = mt_subseq[left:right]

                    # print(f'epilen: {epilen}')
                    # print(f'left: {left}')
                    # print(f'right: {right}')
                    # print(f'wt_subseq_adj: {wt_subseq_adj}')
                    # print(f'mt_subseq_adj: {mt_subseq_adj}')

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

                subseqs.append(entries)

            # its necessary to close the files
            fh.close()
            for epilen in epilens:
                fh_wt[epilen].close()
                fh_mt[epilen].close()
        
            wt_affinities = self.collect_binding_affinities(mhcI_alleles, wt_fname, epilens, 'wt', threads)
            mt_affinities = self.collect_binding_affinities(mhcI_alleles, mt_fname, epilens, 'mt', threads)

            outfile = open(os.path.join(outdir, 'binding_affinities.tsv'),'w')
            BindingAffinities.print_header(outfile)

            for entry in subseqs:
                final_result = {}
                final_result['chrom'] = entry[0]
                final_result['start'] = entry[1]
                final_result['end'] = entry[2]
                final_result['gene_name'] = entry[4]
                final_result['gene_id'] = entry[3]
                final_result['transcript_id'] = entry[5]
                final_result['source'] = entry[9]
                final_result['group'] = entry[10]
                final_result['variant_type'] = entry[11]
                final_result['wt_subseq'] = entry[12]
                final_result['mt_subseq'] = entry[13]
                final_result['new_var_start_pos'] = int(entry[14])
                final_result['vaf'] = float(entry[17])
                final_result['mt_allele_depth'] = int(entry[18])
                final_result['read_depth'] = int(entry[19])
            
                # check results of different epilens
                for epilen_idx in range(0,len(epilens)):

#                    print(entry)

                    # retrieve sequence numbers (from binding affinities results)
                    wt_seqnum = int(entry[24:][epilen_idx*2])
                    mt_seqnum = int(entry[24:][epilen_idx*2+1])

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
                        score = BindingAffinities.calc_ranking_score(final_result['vaf'], final_result['wt_epitope_ic50'], final_result['mt_epitope_ic50'])
                        final_result['ranking_score'] = score
                    
                        BindingAffinities.print_row(final_result, outfile)


            outfile.close()



    @staticmethod
    def collect_binding_affinities(alleles, fnames, epilens, group, threads):
        print("collect binding affinities")
        affinities_results = {}
        no_threads = int(threads)

        # only iterate through keys (alleles)
        # TODO: also incorporate sources
        with concurrent.futures.ThreadPoolExecutor(max_workers=no_threads) as executor:
            futures = {}
            for allele in alleles:
                for epilen in epilens:
                    future = executor.submit(BindingAffinities.calc_binding_affinities, fnames[epilen], allele, epilen, group)
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
    def calc_binding_affinities(fa_file, allele, epilen, group):
        binding_affinities = {}
#    binding_affinities[epilen] = {}
        print(f'calculate binding affinities - allele:{allele} epilen:{epilen} group: {group}')
        call = ['python', 
            'workflow/scripts/mhc_i/src/predict_binding.py', 
            'netmhcpan', 
            allele, 
            str(epilen),
            fa_file]

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
        if wt_ic50 != 'NA':
            fold_change = float(wt_ic50) / float(mt_ic50)
        else:
            fold_change = float(100000.0) / float(mt_ic50)

        score = (1/float(mt_ic50)) + fold_change + vaf * 100.0

        return score

    @staticmethod
    def print_header(outputfile):
        outputfile.write(f"chrom\tstart\tend\tgene_name\tgene_id")
        outputfile.write(f"\ttranscript_id\tsource\tgroup\tvariant_type")
        outputfile.write(f"\tallele\twt_epitope_seq\twt_epitope_ic50\twt_epitope_rank")
        outputfile.write(f"\tmt_epitope_seq\tmt_epitope_ic50\tmt_epitope_rank\tmutation_position")
        outputfile.write(f"\tvaf\tmt_allele_depth\tread_depth\t")
        outputfile.write(f"\tranking_score\n")

    @staticmethod
    def print_row(row, outputfile):
        outputfile.write(f'{row["chrom"]}\t')
        outputfile.write(f'{row["start"]}\t')
        outputfile.write(f'{row["end"]}\t')
        outputfile.write(f'{row["gene_name"]}\t')
        outputfile.write(f'{row["gene_id"]}\t')
        outputfile.write(f'{row["transcript_id"]}\t')
        outputfile.write(f'{row["source"]}\t')
        outputfile.write(f'{row["group"]}\t')
        outputfile.write(f'{row["variant_type"]}\t')
        outputfile.write(f'{row["allele"]}\t')
        outputfile.write(f'{row["wt_epitope_seq"]}\t')
        outputfile.write(f'{row["wt_epitope_ic50"]}\t')
        outputfile.write(f'{row["wt_epitope_rank"]}\t')
        outputfile.write(f'{row["mt_epitope_seq"]}\t')
        outputfile.write(f'{row["mt_epitope_ic50"]}\t')
        outputfile.write(f'{row["mt_epitope_rank"]}\t')
        outputfile.write(f'{row["mutation_position"]}\t')
        outputfile.write(f'{row["vaf"]}\t')
        outputfile.write(f'{row["mt_allele_depth"]}\t')
        outputfile.write(f'{row["read_depth"]}\t')
        outputfile.write(f'{row["ranking_score"]}\n')
        
