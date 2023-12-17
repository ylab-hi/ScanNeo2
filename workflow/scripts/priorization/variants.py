from pathlib import Path
import re
import vcfpy
import sys


class VariantEffectsEntry:
    # create an empty variant effects entry
    def __init__(self):
        self.data = {}
        self.data['chrom'] = ""
        self.data['start'] = ""
        self.data['end'] = ""
        self.data['gene_id'] = ""
        self.data['gene_name'] = ""
        self.data['transcript_id'] = ""
        self.data['transcript'] = ""
        self.data['cds'] = ""
        self.data['cds_breakpoint'] = ""
        self.data['source'] = ""
        self.data['group'] = ""
        self.data['var_type'] = ""
        self.data['wt_subseq'] = ""
        self.data['mt_subseq'] = ""
        self.data["var_start"] = ""
        self.data["local_var_start"] = ""
        self.data["local_var_end"] = ""
        self.data['vaf'] = ""
        self.data['ao'] = ""
        self.data['dp'] = ""
        self.data['NMD'] = "."
        self.data['PTC_dist_ejc'] = "."
        self.data['PTC_exon_number'] = "."
        self.data['NMD_esc_rule'] = "."

    def __init__(self, chrom, start, end, gene_id, gene_name,
                 transcript_id, transcript, cds, cds_breakpoint, source, 
                 group, var_type, wt_subseq, mt_subseq, var_start,
                 local_var_start, local_var_end, vaf, ao, dp):
        self.data = {}
        self.data['chrom'] = chrom
        self.data['start'] = start
        self.data['end'] = end
        self.data['gene_id'] = gene_id
        self.data['gene_name'] = gene_name
        self.data['transcript_id'] = transcript_id
        self.data['transcript'] = transcript
        self.data['cds'] = cds
        self.data['cds_breakpoint'] = cds_breakpoint
        self.data['source'] = source
        self.data['group'] = group
        self.data['var_type'] = var_type
        self.data['wt_subseq'] = wt_subseq
        self.data['mt_subseq'] = mt_subseq
        self.data["var_start"] = var_start
        self.data["local_var_start"] = local_var_start
        self.data["local_var_end"] = local_var_end
        self.data['vaf'] = vaf
        self.data['ao'] = ao
        self.data['dp'] = dp
        self.data['NMD'] = "."
        self.data['PTC_dist_ejc'] = "."
        self.data['PTC_exon_number'] = "."
        self.data['NMD_esc_rule'] = "."


    # getter & setter
    def get_nmd(self):
        return self.data['NMD']

    def set_nmd(self, nmd):
        self.data['NMD'] = nmd


    def scan_nmd(self, anno):
        # check for NMD in framshift events
        if "frameshift" in self.data["var_type"]:

            if self.data["transcript"] == ".":
                self.data["NMD"] = "unknown"
                return 
            else:
                seg2_start = self.data["start"].split('|')[1]
                # search for premature stop codon in cds
                stop_pos, stop_coord = self.find_stop_codon(self.data["cds"],
                                                            self.data["cds_breakpoint"],
                                                            seg2_start)

                # check if stop codon is PTC
                if stop_pos != -1:
                    tid2 = self.data["transcript_id"].split('|')[1]

                    if tid2 == '.':
                        self.data['NMD'] = "unknown"
                        return

                    exoninfo_tid2 = anno[tid2]
                    exons_tid2 = list(exoninfo_tid2.keys())

                    exon_num, dist_ejc = self.annotate_stop_codon(exoninfo_tid2, stop_coord)
                    if exon_num != -1:
                        self.data["PTC_exon_number"] = f'{exon_num}'
                        if max(exons_tid2) > 1:
                            self.data["PTC_dist_ejc"] = dist_ejc

                        # check if NMD is escaped
                        nmd_escape = self.check_escape(exoninfo_tid2,
                                                       stop_coord,
                                                       exon_num,
                                                       dist_ejc)

                        if nmd_escape != -1:
                            self.data["NMD"] = "NMD_escaping_variant"
                            self.data["NMD_esc_rule"] = nmd_escape
                        else:
                            self.data["NMD"] = "NMD_variant"



    def find_stop_codon(self, cds, breakpoint, seg2_start):
        """search for stop codon in coding sequence
        returns the index of the stop codon and 
        genomic coordinate of the end of the second segment
        """
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


    def annotate_stop_codon(self, exoninfo, stop_coord):
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


    def check_escape(self, exoninfo, stop_coord, exon_num, dist_ejc):
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


    def print_vee(self):
        for k,v in self.data.items():
            print(f'{k}: {v}\t', end='')
        print('\n')




class Fusions:
    def __init__(self, fusionsFile, confidence, proteome, annotation, effects):
#        self.variantEffects = VariantEffectsEntry()

        group = Path(fusionsFile).stem.split('_fusions')[0]
        confidenceLevel = {"low": ["low", "medium", "high"],
                           "medium": ["medium", "high"],
                           "high": ["high"]}


        with open(fusionsFile) as fh:
            next(fh) # ignore header
            for line in fh:
                cols = line.rstrip().split('\t')
                if cols[14] in confidenceLevel[confidence]:

                    eventType = cols[8]
                    readingFrame = cols[15]
                    consequence = self.determineConsequence(readingFrame, eventType)
                    if consequence == None:
                        continue

                    peptideSeq = cols[28]
                    if peptideSeq == '.':
                        continue
                    
                    # get genomic coordiantes
                    chrom, start, end = self.determineGenomicCoordinates(cols[4],cols[5])
                    gene_name = cols[0] + '|' + cols[1]
                    gene_id = cols[20] + '|' + cols[21]

                    tid1 = cols[22]
                    tid2 = cols[23]

                    transcript_id = tid1 + '|' + tid2

                    # print(f'tid2: {tid2}')
                    # print(f'transcript_id: {transcript_id}')

                    ao = str(len(cols[29]))
                    dp = str(cols[12]+cols[13])

                    transcript = cols[27].upper()
                    cds, cds_breakpoint = self.determineCodingSequence(transcript)

                    # retrieve polypeptide sequence of 
                    fusion_pepseq = cols[28]
                    if fusion_pepseq == '.':
                        continue

                    # determine (complete) wildtype peptide sequence (which is the wt of first gene)
                    gene1_wt_seq = '.'
                    # retrieve 
                    if tid1 in proteome:
                        gene1_wt_seq = proteome[tid1]
                    wt_subseq, mt_subseq, var_start, local_var_start, local_var_end = self.determine_peptide_sequence(fusion_pepseq, 
                                                                                                                      gene1_wt_seq)
                    vee = VariantEffectsEntry(chrom, start, end, gene_id,
                                              gene_name, transcript_id,
                                              transcript, cds, cds_breakpoint, 
                                              "fusion", group, consequence, # var_type
                                              wt_subseq, mt_subseq, var_start, 
                                              local_var_start, local_var_end,
                                              ".", # vaf not available for fusions 
                                              ao, dp)

                    #vee.print_vee()
                    vee.scan_nmd(annotation)
                    effects.print_entry(vee)


    def determineCodingSequence(self, transcript):
        """ determines the coding sequence """
#        print(f'transcript: {transcript}')

        #search for start codon in transcript
        start_codon = re.search(r'ATG', transcript)
#        print(f'start_codon: {start_codon}')
        if start_codon: # start codon could have been found
#            print(f'start_codon: {start_codon}')
            cds = transcript[start_codon.start()+3:] # determine the coding sequence
            cds_bp = cds.find('|')
            cds = cds.replace('|', '')
            return cds, cds_bp

        else:
            return '.', '.'


    def determineConsequence(self, readingFrame, event):
        consequence = ""
        if readingFrame == '.' or readingFrame == "stop_codon":
            return None
        else:
            if readingFrame == "in-frame":
                consequence = "inframe_"
            elif readingFrame == "out-of-frame":
                consequence = "frameshift_"

            evts = event.split('/')
            for evt in evts:
                if (evt == "read-through" or
                    evt == "ITD" or
                    evt == "5'-5'" or
                    evt == "3'-3'"):
                    return None
                elif evt == "translocation":
                    consequence += "trs"
                    return consequence
                elif evt == "duplication":
                    consequence += "dup"
                    return consequence
                elif evt == "deletion":
                    consequence += "del"
                    return consequence
                elif evt == "inversion":
                    consequence += "inv"
                    return consequence
        return None


    def determineGenomicCoordinates(self, breakpoint1, breakpoint2):
        chrom1 = breakpoint1.split(':')[0]
        chrom2 = breakpoint2.split(':')[0]

        start1 = breakpoint1.split(':')[1]
        start2 = breakpoint2.split(':')[1]

        chrom = chrom1 + '|' + chrom2
        start = start1 + '|' + start2  
        end = -1

        return chrom, start, end

    def determine_peptide_sequence(self, fusion_seq, wildtype_seq):
        fusion_seg1 = fusion_seq.split('|')[0]
        fusion_seg2 = fusion_seq.split('|')[1]
        var_start = len(fusion_seg1) + 1 # position of the fusion

        # extract sequences
        fusion_seg1_sub = fusion_seg1.split('?')[-1]
        fusion_seg2_sub = fusion_seg2.split('?')[0].split('*')[0]

        # mt sequence
        if len(fusion_seg1_sub) >= 24:
            mt = fusion_seg1_sub[-24:] 
            local_var_start = 25 # first base of fused transcript
            if len(fusion_seg2_sub) >= 24:
                mt += fusion_seg2_sub[:24]
            else:
                mt += fusion_seg2_sub
        elif len(fusion_seg1_sub) < 24:
            mt = fusion_seg1_sub
            local_var_start = len(fusion_seg1_sub) + 1
            if len(fusion_seg2_sub) >= 24:
                mt += fusion_seg2_sub[:24]
            else:
                mt += fusion_seg2_sub
       
        # arriba returns variants (in lower case)
        mt = mt.upper()

        # end of variant is always the length of the sequence
        local_var_end = len(mt)

        # wt sequences
        wt = '.'
        # assume that the fusion sequence is long 

        if len(fusion_seg1_sub) > 10:
            # search for variants in the first segment 
            variants = find_all_lowercase(fusion_seg1_sub)
            if variants != []: # variants within the segment
                # extract sequence that is wildtype (in fusion - before first lowercase)
                wt_prefix = fusion_seg1_sub[:variants[0]]
            else:
                wt_prefix = fusion_seg1_sub

            wt = self.extract_wildtype_sequence(wildtype_seq, wt_prefix, len(mt))
        return wt, mt, var_start, local_var_start, local_var_end


    # extract the wildtype sequence (using ensemble info and fusion sequence)
    def extract_wildtype_sequence(self,wt_seq, wt_prefix, mt_len):
        # search for subsequence in wildtype peptide sequences
        wt_start = wt_seq.find(wt_prefix)
        if wt_start != -1: # wildtype sequence (before first variant) can be found in ensemble
            return self.fillup_wildtype_sequence(wt_seq[wt_start:], mt_len)
        else: # wildrtype sequence cannot be found in ensemble
            if wt_prefix.isupper():
                return self.fillup_wildtype_sequence(wt_prefix, mt_len)
            else: # if not all in uppercase wt cannot be determined
                return '.'


    # makes sure that length of wildtype sequence equals length of mutant sequence
    def fillup_wildtype_sequence(self, wt_seq,mt_len):
        # wildtype sequence of seg1 exceeds length of mt (take so much until length is reached)
        if len(wt_seq) >= mt_len:
            wt = wt_seq[:mt_len]
        else:
            wt = wt_seq + '$'*(mt_len-len(wt_seq))
        return wt



class Variants():
    def __init__(self, variants, effects):

        input_vcfs = variants.split(' ')
        for vcf in input_vcfs:
            vcf_reader = vcfpy.Reader(open(vcf, 'r'))
            csq_format = self.parse_csq_format(vcf_reader.header)
        
            transcript_count = {}
            for entry in vcf_reader:
                # FILTER (when applicable)
                if (entry.INFO['SRC'] == 'snv' or
                    entry.INFO['SRC'] == 'short_indel'):
                    if 'PASS' not in entry.FILTER:
                        continue
            
                # resolved alleles specific to vep
                alleles_vep = self.resolve_alleles(entry)

                chromosome = entry.CHROM
                start = entry.affected_start
                stop = entry.affected_end
                ref = entry.REF
                alts = entry.ALT
            
                for alt_i, alt_v in enumerate(alts):
                # for alt in alts:
                    csq_allele = alleles_vep[str(alt_v.value)]
                    csq_fields = self.parse_csq_entries(
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

                    nmd = '.'
                    csq = self.resolve_consequence(field['Consequence'])
                    if csq is None:
                        continue
                    elif csq == "frameshift":
                        if field["NMD"] != 'NMD_escaping_variant':
                            continue
                        elif field["NMD"] == "NMD_escaping_variant":
                            nmd = "NMD_escaping_variant"

                        if field["DownstreamProtein"] == "":
                            continue
                    
                    vaf = self.determine_variant_allele_frequency(entry, alt_i)
                    # determine allele depth (number of reads)
                    mt_ad = self.determine_allele_depth(entry, alt_i)
                    dp = self.determine_read_depth(entry)
                    
                    if field["Amino_acids"]:
                        aa_change = field["Amino_acids"]
                    else:
                        continue

                    if csq == "frameshift":
                        var_start = self.get_variant_startpos(field['Protein_position'])
                        wt_seq = field["WildtypeProtein"]
                        # mutant/variant sequence is wt until mutation start
                        mt_seq = wt_seq[:var_start] + field["DownstreamProtein"]
                        # length of variant is length of downstream sequence
                        var_len = len(field["WildtypeProtein"])
                        
                        # makes sure that wt and mt are of same length
                        wt_subseq, mt_subseq, n_var_start = self.determine_fs_subseq(
                                wt_seq,
                                mt_seq,
                                var_start
                        )
                        
                        # print('-------frameshift')
                        # print(f'field-Protein_Position:{field["Protein_position"]}')
                        # print(f'Downstream_Protein: {field["DownstreamProtein"]}')
                        # print(f'var_start: {var_start}')
                        # print(f'aa_change:{aa_change}')
                        # print(f'wt_seq: {wt_seq}')
                        # print(f'mt_seq: {mt_seq}')

                        # print(f'wt_subseq: {wt_subseq}')
                        # print(f'mt_subseq: {mt_subseq}')
                        # print(f'n_var_start: {n_var_start}')

                    
                        vee = VariantEffectsEntry(entry.CHROM, 
                                                  start, 
                                                  stop,
                                                  gene_id,
                                                  gene_name,
                                                  transcript_id,
                                                  '.', # transcript
                                                  '.', # cds
                                                  '.', # cds_breakpoint
                                                  entry.INFO["SRC"],
                                                  entry.INFO["GRP"],
                                                  csq,
                                                  wt_subseq,
                                                  mt_subseq,
                                                  var_start,
                                                  n_var_start,
                                                  len(mt_subseq)-1,
                                                  vaf,
                                                  mt_ad,
                                                  dp)

                        vee.set_nmd(nmd)
                        effects.print_entry(vee)



                    elif (csq == "missense" or
                          csq == "inframe_ins" or 
                          csq == "inframe_del"):
                        

                        # retrieve the start and end of the variant / initial variant start 
                        var_start = self.get_variant_startpos(field["Protein_position"])
                        wt_seq = field["WildtypeProtein"] # wildtype peptide sequence

                        wt_aa_change, mt_aa_change = self.determine_aa_change(aa_change)
                        
                        # scan for stop codons
                        wt_aa_change, wt_stop_codon = self.scan_stop_codon(wt_aa_change)
                        mt_aa_change, mt_stop_codon = self.scan_stop_codon(mt_aa_change)

                        mt_seq = wt_seq[:var_start] + mt_aa_change
                        if not mt_stop_codon:
                            mt_seq += wt_seq[var_start+len(wt_aa_change):]

                        # change end when deletion and insertion (as opposed to missense)
#                            var_end = var_start
                        # if 'del' in csq:
                        # elif 'ins' in csq:
                            # var_end = var_start + len(mt_aa_change)
                        # else:
                            # var_end = var_start

                        # print(f'var_start: {var_start}')
                        # print(f'var_end: {var_end}')
                        # print('-------')
                        # print(f'field-Protein_Position:{field["Protein_position"]}')
                        # print(f'aa_change:{aa_change}')
                        # print(f'wt_seq: {wt_seq}')

                        # print(f'wt_aa_change: {wt_aa_change}')
                        # print(f'mt_aa_change: {mt_aa_change}')

                        # print(f'mt_seq: {mt_seq}')

                        # # generate subsequence
                        wt_subseq, mt_subseq, n_var_start, n_var_end = self.determine_subseq(wt_seq,
                                                                     mt_seq, 
                                                                     var_start,
                                                                     len(mt_aa_change))


                        vee = VariantEffectsEntry(entry.CHROM, 
                                                  start, 
                                                  stop,
                                                  gene_id,
                                                  gene_name,
                                                  transcript_id,
                                                  '.', # transcript
                                                  '.', # cds
                                                  '.', # cds_breakpoint
                                                  entry.INFO["SRC"],
                                                  entry.INFO["GRP"],
                                                  csq,
                                                  wt_subseq,
                                                  mt_subseq,
                                                  var_start,
                                                  n_var_start,
                                                  n_var_end,
                                                  vaf,
                                                  mt_ad,
                                                  dp)

                        effects.print_entry(vee)




                        # print(f'wt_subseq: {wt_subseq}')
                        # print(f'mt_subseq: {mt_subseq}')
                        # print(f'n_var_start: {n_var_start}')
                        # print(f'n_var_end: {n_var_end}')



                        #output_row["aa_change"] = field["Amino_acids"]
                    # else:
                        # output_row["aa_change"] = "NA"
                        # continue

    @staticmethod
    def determine_subseq(wt_seq, mt_seq, var_start, var_len):
        """ determine the subsequence to consider for priorization """
        # left and right are start/end of the subsequence
        # new_var_start is the new start of the variant (after subsetting) 
        if var_start <= 24 :
            left = 0 # start of subsequence
            new_var_start = var_start
        else:
            left = var_start - 24
            new_var_start = 24

        if var_start + var_len + 24 >= len(mt_seq):
            right = len(mt_seq)
        else:
            right = var_start  + var_len + 24
        new_var_end = new_var_start + var_len

        # get wildtype sequence with new subsets
        mt_subseq = mt_seq[left:right] 
        if right >= len(wt_seq):
            wt_subseq = wt_seq[left:]
        else:
            wt_subseq = wt_seq[left:right]

        if wt_subseq < mt_subseq:
            for i in range(len(wt_subseq),len(mt_subseq)):
                wt_subseq += '$'

        return wt_subseq, mt_subseq, new_var_start, new_var_end

    @staticmethod
    def determine_fs_subseq(wt_seq, mt_seq, start_pos_var):
        subseqs = {}
        if start_pos_var < 24:
            subseq_start = 0
        else:
            subseq_start = start_pos_var - 24

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


    @staticmethod
    def scan_stop_codon(sequence):
        """ checks for stop codons and removes them """
        stop_codon = False
        if '*' in sequence:
            sequence = sequence.split('*')[0]
            stop_codon = True
        if 'X' in sequence:
            sequence = sequence.split('X')[0]
            stop_codon = True

        return sequence, stop_codon


    @staticmethod
    def determine_aa_change(aa_change):
        wt_aa_change, mt_aa_change = aa_change.split('/')
        if wt_aa_change == '-':
            wt_aa_change = ''
        if mt_aa_change == '-':
            mt_aa_change = ''

        return wt_aa_change, mt_aa_change



    def get_variant_startpos(self, protein_position):
        if '-' not in protein_position: # single number
            startpos = int(protein_position) - 1
        else:
            if protein_position.split('/')[0] == '-':
                startpos = int(protein_position.split('-')[0])
            else:
                startpos = int(protein_position.split('-')[0]) - 1

        return startpos



    def parse_csq_entries(self, csq_entries, csq_format, csq_allele):
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


    def parse_csq_format(self,vcf_header):
        if vcf_header.get_info_field_info('CSQ') is None:
            sys.exit("Failed to extract format string form info description for tag (CSQ)")
        else:
            csq_header = vcf_header.get_info_field_info('CSQ')
            format_pattern = re.compile("Format: (.*)")
            match = format_pattern.search(csq_header.description)
            return match.group(1)


    def resolve_alleles(self, entry):
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


    def resolve_consequence(self, consequence_string):
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

    def determine_variant_allele_frequency(self, entry, i):
        vaf = -1
        if (entry.INFO['SRC'] == 'short_indel' or 
            entry.INFO['SRC'] == 'snv'):
            vaf = entry.calls[0].data['AF'][i]

        if entry.INFO['SRC'] == 'long_indel':
            vaf = entry.INFO['AB']

        return vaf

    def determine_read_depth(self, entry):
        dp = -1
        if (entry.INFO['SRC'] == 'short_indel' or 
            entry.INFO['SRC'] == 'snv'):
            dp = entry.calls[0].data['DP']

        if (entry.INFO['SRC'] == 'long_indel' or 
            entry.INFO['SRC'] == 'exitron'):
            dp = entry.INFO['DP']

        return dp


    def determine_allele_depth(self, entry, i):
        #wt_ad = -1
        mt_ad = -1

        # GATK
        if (entry.INFO['SRC'] == 'short_indel' or 
            entry.INFO['SRC'] == 'snv'):
            #wt_ad = entry.calls[0].data['AD'][0]
            mt_ad = entry.calls[0].data['AD'][i+1]

        # transIndel
        if (entry.INFO['SRC'] == 'long_indel' or 
            entry.INFO['SRC'] == 'exitron'):
            mt_ad = entry.INFO['AO']
            #wt_ad = str(int(float(entry.INFO['DP']) - float(mt_ad)))

        return mt_ad


# augmented transcript
class VariantEffects:
    def __init__(self, output_dir):

        output_dir_path =  Path(output_dir)
        if not output_dir_path.exists():
            output_dir_path.mkdir()

        self.variantEffectsFile = Path(output_dir_path, "variant_effects.tsv")
        # create dir if it not exists

        self.fh = open(self.variantEffectsFile, 'w')



    # TODO: no hard coding of output
    def print_header(self):
        self.fh.write("chrom\tstart\tend\tgene_id\tgene_name\ttranscript_id\t")
        self.fh.write("transcript\tcds\tcds_breakpoint\tsource\tgroup\t")
        self.fh.write("var_type\twt_subseq\tmt_subseq\tvar_start\t")
        self.fh.write("local_var_start\tlocal_var_end\t")
        self.fh.write("vaf\tao\tdp\tNMD\tPTC_dist_ejc\tPTC_exon_number\t")
        self.fh.write("NMD_esc_rule\n")

    def print_entry(self, entry):
        self.fh.write(f'{entry.data["chrom"]}\t')
        self.fh.write(f'{entry.data["start"]}\t')
        self.fh.write(f'{entry.data["end"]}\t')
        self.fh.write(f'{entry.data["gene_id"]}\t')
        self.fh.write(f'{entry.data["gene_name"]}\t')
        self.fh.write(f'{entry.data["transcript_id"]}\t')
        self.fh.write(f'{entry.data["transcript"]}\t')
        self.fh.write(f'{entry.data["cds"]}\t')
        self.fh.write(f'{entry.data["cds_breakpoint"]}\t')
        self.fh.write(f'{entry.data["source"]}\t')
        self.fh.write(f'{entry.data["group"]}\t')
        self.fh.write(f'{entry.data["var_type"]}\t')
        self.fh.write(f'{entry.data["wt_subseq"]}\t')
        self.fh.write(f'{entry.data["mt_subseq"]}\t')
        self.fh.write(f'{entry.data["var_start"]}\t')
        self.fh.write(f'{entry.data["local_var_start"]}\t')
        self.fh.write(f'{entry.data["local_var_end"]}\t')
        self.fh.write(f'{entry.data["vaf"]}\t')
        self.fh.write(f'{entry.data["ao"]}\t')
        self.fh.write(f'{entry.data["dp"]}\t')
        self.fh.write(f'{entry.data["NMD"]}\t')
        self.fh.write(f'{entry.data["PTC_dist_ejc"]}\t')
        self.fh.write(f'{entry.data["PTC_exon_number"]}\t')
        self.fh.write(f'{entry.data["NMD_esc_rule"]}\n')
        

    def close_file(self):
        self.fh.close()


# deter
def find_all_lowercase(seq):
    matches = []
    for i in range(len(seq)):
        if seq[i].islower():
            matches.append(i)
    return matches
