"""MHC binding-affinity prediction.

Extracts candidate epitope subsequences from the variant-effects table, runs netMHCpan / netMHCIIpan
(IEDB tools) in batches across the sample's HLA alleles and the configured epitope lengths, and writes
the per-variant-type neoepitope tables with binding affinity, rank, and agretopicity.
"""

import tempfile
import os
import contextlib
import concurrent.futures
import subprocess
from pathlib import Path
import utility as ut

BATCH_SIZE = 500
PREDICTION_TIMEOUT_SEC = 3600  # per-batch wall-clock cap for netMHCpan / netMHCIIpan

class BindingAffinities:
    def __init__(self, threads):
        self.threads = threads
        pass

    def start(self, allele_file, epitope_lengths, output_dir, mhc_class, vartype):
        # create temorary_directory
        with tempfile.TemporaryDirectory() as tmp_seqs:
            self.get_alleles(allele_file)

            # initialize filenames
            wt_fname = {}
            mt_fname = {}

            # counters to store number of sequence in fa
            wt_cnt = {}
            mt_cnt = {}

            epilens = self.extract_epilens(epitope_lengths)

            for epilen in epilens:
                wt_fname[epilen] = os.path.join(tmp_seqs, f'wt_{epilen}.fa')
                mt_fname[epilen] = os.path.join(tmp_seqs, f'mt_{epilen}.fa')

                # initialise counter
                wt_cnt[epilen] = 1
                mt_cnt[epilen] = 1

            subseqs = []

            # Open all per-length write handles plus the variant effects input
            # together; ExitStack closes (and flushes) them when the block
            # exits, before the binding affinity subprocess reads them.
            with contextlib.ExitStack() as stack:
                fh_wt = {epilen: stack.enter_context(open(wt_fname[epilen], 'w')) for epilen in epilens}
                fh_mt = {epilen: stack.enter_context(open(mt_fname[epilen], 'w')) for epilen in epilens}
                fh = stack.enter_context(open(Path(output_dir, f"{vartype}_variant_effects.tsv"), 'r'))
                next(fh)   # skip header
                for line in fh:
                    entries = line.rstrip().split('\t')
                    # drop the wt-padding sentinel: '$' marks positions where
                    # the wildtype has no residue (mt is longer). It is internal
                    # to effects.py and must not reach epitope output or scoring.
                    entries[9] = entries[9].split('$')[0]
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

                        wt_subseq_adj = wt_subseq[left:right+1]
                        mt_subseq_adj = mt_subseq[left:right+1]

                        # determine the epitope sequences ('$' was already
                        # stripped from wt_subseq when the tsv line was parsed)
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

            total_seqs = max((wt_cnt.get(epilens[0], 1),
                              mt_cnt.get(epilens[0], 1))) - 1
            print(f"calculate binding affinities for {total_seqs} sequences "
                  f"({len(self.alleles)} alleles, epitope lengths: "
                  f"{','.join(map(str, epilens))})...", flush=True)

            wt_affinities, mt_affinities = self.collect_binding_affinities(
                self.alleles,
                {'wt': wt_fname, 'mt': mt_fname},
                epilens,
                mhc_class,
                self.threads)
            print("Done", flush=True)
            
            with open(os.path.join(output_dir,
                                   f"{vartype}_{mhc_class}_neoepitopes.txt"), "w") as outfile:
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
                            # find every position of the epitope within
                            # mt_subseq -- the same k-mer can occur more than
                            # once, and find() alone returns only the first
                            positions = []
                            pos = mt_subseq.find(epitope)
                            while pos != -1:
                                positions.append(pos)
                                pos = mt_subseq.find(epitope, pos + 1)

                            # keep the epitope only if some occurrence spans the
                            # variant [aa_var_start, aa_var_end] -- one lying
                            # entirely up-/downstream is pure wildtype, not a
                            # neoepitope. Use that occurrence for the wt epitope.
                            startpos = next((p for p in positions
                                             if p <= aa_var_end
                                             and p + len(epitope) - 1 >= aa_var_start),
                                            None)
                            if startpos is None:
                                continue

                            # mt[epitope] = (allele, start, end, ic50, rank)
                            final["mt_epitope_seq"] = epitope
                            final["allele"] = mt[epitope][0]
                            final["mt_epitope_ic50"] = mt[epitope][3]
                            final["mt_epitope_rank"] = mt[epitope][4]

                            # the wt epitope occupies the same coordinates in
                            # wt_subseq as the mt epitope does in mt_subseq
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



    def get_alleles(self, allele_file):
        self.alleles = {}
        with open(allele_file, 'r') as fh_alleles:
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
    def collect_binding_affinities(alleles, fnames, epilens, mhc_class, threads):
        """Run binding-affinity prediction for wt and mt over every allele,
        epitope length and FASTA batch in a single thread pool.

        fnames is {'wt': {epilen: path}, 'mt': {epilen: path}}. Returns
        (wt_affinities, mt_affinities), each
        {epilen: {global_seqnum: {epitope: tuple}}}.
        """
        affinities = {grp: {epilen: {} for epilen in epilens}
                      for grp in ('wt', 'mt')}

        with contextlib.ExitStack() as stack:
            # split every (group, epilen) FASTA up front and enumerate the
            # individual (group, allele, epilen, batch) work units
            units = []
            for group in ('wt', 'mt'):
                for epilen in epilens:
                    fa_file = fnames[group][epilen]
                    if os.stat(fa_file).st_size == 0:
                        continue
                    batch_dir = stack.enter_context(tempfile.TemporaryDirectory())
                    batches = BindingAffinities.split_fasta_into_batches(
                        fa_file, batch_dir)
                    for allele in alleles:
                        for batch_file, offset in batches:
                            units.append(
                                (group, allele, epilen, batch_file, offset))

            print(f"  submitting {len(units)} prediction jobs "
                  f"({len(alleles)} alleles x {len(epilens)} epitope lengths "
                  f"x wt/mt x FASTA batches)", flush=True)

            # one pool over all units -- the prediction tool numbers each batch
            # file from 1, so offset translates that to a global seqnum
            completed = 0
            with concurrent.futures.ThreadPoolExecutor(
                    max_workers=int(threads)) as executor:
                futures = {}
                for group, allele, epilen, batch_file, offset in units:
                    call = BindingAffinities._build_call(
                        batch_file, allele, epilen, mhc_class)
                    future = executor.submit(BindingAffinities._run_prediction,
                                             call, batch_file, group, mhc_class)
                    futures[future] = (group, epilen, offset)

                for future in concurrent.futures.as_completed(futures):
                    group, epilen, offset = futures[future]
                    completed += 1
                    print(f"  [{completed}/{len(units)}] completed", flush=True)
                    dest = affinities[group][epilen]
                    for seqnum, epitopes in future.result().items():
                        global_seqnum = offset + seqnum
                        if global_seqnum not in dest:
                            dest[global_seqnum] = epitopes
                        else:
                            for seq, val in epitopes.items():
                                dest[global_seqnum].setdefault(seq, val)

        return affinities['wt'], affinities['mt']
    
    @staticmethod
    def split_fasta_into_batches(fa_file, batch_dir, batch_size=BATCH_SIZE):
        """Split a FASTA file into batch files of at most batch_size sequences.

        Returns a list of (batch_file, offset) tuples. The prediction tool
        numbers the sequences in each file it receives starting from 1, so the
        offset -- the count of sequences in all preceding batches -- is needed
        to recover a global sequence number: global = offset + per-batch number.
        """
        batches = []
        current_lines = []
        seq_count = 0       # sequences in the batch currently being built
        total_count = 0     # sequences already flushed to preceding batches

        with open(fa_file, 'r') as fh:
            for line in fh:
                if line.startswith('>'):
                    # start of a new sequence – flush if batch is full
                    if seq_count >= batch_size and current_lines:
                        batch_file = os.path.join(
                            batch_dir, f'batch_{len(batches)}.fa')
                        with open(batch_file, 'w') as bf:
                            bf.writelines(current_lines)
                        batches.append((batch_file, total_count))
                        total_count += seq_count
                        current_lines = []
                        seq_count = 0
                    seq_count += 1
                current_lines.append(line)

        # write remaining sequences
        if current_lines:
            batch_file = os.path.join(
                batch_dir, f'batch_{len(batches)}.fa')
            with open(batch_file, 'w') as bf:
                bf.writelines(current_lines)
            batches.append((batch_file, total_count))

        return batches

    @staticmethod
    def _run_prediction(call, fa_file, group, mhc_class):
        """Run a single prediction subprocess and parse its output into a
        binding_affinities dict keyed by (seqnum -> epitope_seq -> tuple).

        A hung or failing batch is skipped (returns an empty dict) so one
        bad batch does not stall the whole run."""
        binding_affinities = {}

        try:
            result = subprocess.run(call,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True,
                                    timeout=PREDICTION_TIMEOUT_SEC,
                                    check=True)
        except subprocess.TimeoutExpired:
            print(f"    WARNING: prediction timed out after "
                  f"{PREDICTION_TIMEOUT_SEC}s for {fa_file}", flush=True)
            return binding_affinities
        except subprocess.CalledProcessError as e:
            print(f"    WARNING: prediction failed for {fa_file}: "
                  f"{e.stderr}", flush=True)
            return binding_affinities

        predictions = result.stdout.rstrip().split('\n')[1:]

        for line in predictions:
            entries = line.split('\t')
            if group == 'mt':
                if float(entries[8]) >= 500:
                    continue

            # start and end in sequence (0-based)
            start = int(entries[2]) - 1
            end = int(entries[3]) - 1

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
    def _build_call(batch_file, allele, epilen, mhc_class):
        """Build the subprocess command list for a prediction tool."""
        if mhc_class == "mhc-I":
            return ["python",
                    "workflow/scripts/mhc_i/src/predict_binding.py",
                    "netmhcpan", allele, str(epilen), batch_file]
        elif mhc_class == "mhc-II":
            return ["python",
                    "workflow/scripts/mhc_ii/mhc_II_binding.py",
                    "netmhciipan_ba", allele, batch_file, str(epilen)]
    
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
        
    
