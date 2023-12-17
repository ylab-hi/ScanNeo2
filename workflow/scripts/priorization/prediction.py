import tempfile
import pdb
import os
import concurrent.futures

class BindingAffinities:
    def __init__(self, effects, mhcI, mhcII, mhcI_len, mhcII_len):
    
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
                print(entries)
                # print(line)

                for epilen in epilens:
                    local_var_start = int(entries[15])
                    local_var_end = int(entries[16])

                    wt_subseq = entries[12]
                    mt_subseq = entries[13]

                    # adjust length of seqs (according to epilen)
                    if local_var_start >= epilen + 1:
                        left = local_var_start - (epilen + 1)
                    else:
                        left = 0

                    if local_var_end + (epilen - 1) <= len(mt_subseq):
                        right = local_var_end + (epilen - 1)
                    else:
                        right = len(mt_subseq) - 1

                    wt_subseq_adj = wt_subseq[left:right]
                    mt_subseq_adj = mt_subseq[left:right]

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
        
            wt_affinities = self.collect_binding_affinities(alleles, wt_fname, epilens, 'wt', options.threads)
            mt_affinities = self.collect_binding_affinities(alleles, mt_fname, epilens, 'mt', options.threads)


    @staticmethod
    def collect_binding_affinities(alleles, fnames, epilens, group, threads):
        affinities_results = {}
        number_of_threads = int(threads)

        with concurrent.futures.ThreadPoolExecutor(max_workers=number_of_threads) as executor:
            futures = {}
            for allele in alleles:
                for epilen in epilens:
                    future = executor.submit(calc_binding_affinities, fnames[epilen], allele, epilen, group)
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




            

