from pathlib import Path
import re
import vcfpy
import sys
import reference

import pdb


class Fusions:
    def __init__(self, fusionsFile, 
                 confidence, 
                 proteome, 
                 annotation, 
                 variant_effects):
#        self.variantEffects = VariantEffectsEntry()

        group = Path(fusionsFile).stem.split('_fusions')[0]
        confidence_level = {"low": ["low", "medium", "high"],
                           "medium": ["medium", "high"],
                           "high": ["high"]}


        with open(fusionsFile) as fh:
            next(fh) # ignore header
            for line in fh:
                cols = line.rstrip().split('\t')
                if cols[14] in confidence_level[confidence]:
                    event_type = cols[8]
                    reading_frame = cols[15]
                    csq = self.determine_consequence(reading_frame, event_type)
                    if csq == None:
                        continue

                    peptide_seq = cols[28]
                    if peptide_seq == '.': # could be determined
                        continue
                    
                    # get genomic coordiantes
                    chrom, start, end = self.determine_coordinates(cols[4],cols[5])
                    gene_name = cols[0] + '|' + cols[1]
                    gene_id = cols[20] + '|' + cols[21]

                    tid1 = cols[22]
                    tid2 = cols[23]

                    transcript_id = tid1 + '|' + tid2

                    ao = str(len(cols[29]))
                    dp = str(cols[12]+cols[13])

                    # get the fusion transcript
                    transcript = cols[27].upper()
                    for symbol in ['[', ']', '(', ')', '_']:
                        transcript = transcript.replace(symbol, '')
                    transcript_bp = transcript.find('|')
                    transcript = transcript.replace('|', '')

                    # retrieve polypeptide sequence of 
                    fusion_pepseq = cols[28]
                    if fusion_pepseq == '.':
                        continue

                    # determine (complete) wildtype peptide sequence (which is the wt of first gene)
                    gene1_wt_seq = ''
                    # retrieve 
                    if tid1 in proteome:
                        gene1_wt_seq = proteome[tid1]
                    wt_seq, mt_seq, var_start = self.determine_peptide_sequence(fusion_pepseq,
                                                                                gene1_wt_seq)

                    variant_effects.change_entry(chrom=chrom,
                                               start=start,
                                               end=None,
                                               gene_id=gene_id,
                                               gene_name=gene_name,
                                               transcript_id=transcript_id,
                                               transcript=transcript,
                                               transcript_bp=transcript_bp,
                                               source="fusion",
                                               group=group,
                                               var_type=csq,
                                               var_start=var_start,
                                               wt_seq=wt_seq,
                                               mt_seq=mt_seq,
                                               vaf=-1,
                                               ao=ao,
                                               dp=dp,
                                               nmd_event=None)

                    # check if variant is self dissimilar
                    if variant_effects.self_dissimilarity():
                        variant_effects.write_entry()


    @staticmethod
    def determine_cds(transcript):
        """ determines the coding sequence of the transcript """
        start_codon = re.search(r'ATG', transcript)
        if start_codon:
            # determines the cds - starting after start codon
            cds = transcript[start_codon.start()+3:]
            cds_bp = cds.find('|')
            cds = cds.replace('|', '')
            return cds, cds_bp
        else:
            return None, None


    @staticmethod
    def determine_consequence(reading_frame, event_type):
        """ determines the consequence of the fusion event """
        consequence = ""
        # stop_codon indicates that 
        if reading_frame == '.' or reading_frame == "stop_codon":
            return None
        else:
            if reading_frame == "in-frame":
                consequence = "inframe_"
            elif reading_frame == "out-of-frame":
                consequence = "frameshift_"

            evts = event_type.split('/')
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

    @staticmethod
    def determine_coordinates(breakpoint1, breakpoint2):
        chrom1 = breakpoint1.split(':')[0]
        chrom2 = breakpoint2.split(':')[0]

        start1 = breakpoint1.split(':')[1]
        start2 = breakpoint2.split(':')[1]

        chrom = chrom1 + '|' + chrom2
        start = start1 + '|' + start2  
        end = -1

        return chrom, start, end

    @staticmethod
    def determine_peptide_sequence(fusion_seq, wildtype_seq):
        # extract sequence - without undefined seq (e.g., ?) 
        fusion_seg1 = fusion_seq.split('|')[0].split('?')[-1]
        fusion_seg2 = fusion_seq.split('|')[1].split('?')[0].split('*')[0]

        # first position of the 2nd fusion transcript is the var_start
        var_start = len(fusion_seg1)

        # arriba returns variants (in lower case)
        mt =  (fusion_seg1 + fusion_seg2).upper()
        """We require the first segment to be at least 10 bases long. The 
        reason is that for MHCI we need at least 7 bases (1 base in seg2)"""
        if len(fusion_seg1) >= 10:
            """search for variants in the first segments (lowercase) - we need
            the wildtype bases to determine the full wildtype sequence"""
            variants = find_all_lowercase(fusion_seg1)
            if variants != []:
                # extract sequence that is wildtype (in fusion - before lowercase)
                wt_prefix = fusion_seg1[:variants[0]]
            else: # all of seg1 is wildtype
                wt_prefix = fusion_seg1

            wt = Fusions.determine_wildtype_sequence(wildtype_seq, wt_prefix)
        else:
            wt = ''

        return wt, mt, var_start


    # extract the wildtype sequence (using ensemble info and fusion sequence)
    @staticmethod
    def determine_wildtype_sequence(wt_seq, wt_prefix):
        # search for subsequence in wildtype peptide sequences
        wt_start = wt_seq.find(wt_prefix)
        if wt_start != -1: # wildtype sequence (before first variant) can be found in ensemble
            return wt_seq[wt_start:]
        else: # wildrtype sequence cannot be found in ensemble and seg1 is incomplete
            return ''


class Variants():
    def __init__(self, variants_input, annotation, variant_effects):

        input_vcfs = variants_input.split(' ')
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

                chrom = entry.CHROM
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

                    nmd = None
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
#                        var_len = len(field["WildtypeProtein"])

                        stop_pos = mt_seq.find('*')
                        unknown_pos = mt_seq.find('X')

                        if stop_pos != -1:
                            mt_seq = mt_seq[:stop_pos]

                        if unknown_pos != -1:
                            mt_seq = mt_seq[:unknown_pos]


                    elif (csq == "missense" or
                          csq == "inframe_INS" or 
                          csq == "inframe_DEL"):
                        
                        # retrieve the start and end of the variant / initial variant start 
                        var_start = self.get_variant_startpos(field["Protein_position"])
                        wt_aa_change, mt_aa_change = self.determine_aa_change(aa_change)
                        
                        # scan for stop codons
                        wt_aa_change, wt_stop_codon = self.scan_stop_codon(wt_aa_change)
                        mt_aa_change, mt_stop_codon = self.scan_stop_codon(mt_aa_change)

                        wt_seq = field["WildtypeProtein"] # wildtype peptide sequence
                        # check if there are unknown amino Amino_acids
                        if 'X' in wt_seq:
                            # needs to occur before varstart...
                            unknown_pos = wt_seq.find('X')
                            if unknown_pos != -1 and unknown_pos < var_start:
                                var_start = var_start - unknown_pos - 1
                                wt_seq = wt_seq[unknown_pos+1:]
                            else:
                                # ...otherwise subsequence is altered - skip
                                continue

                        mt_seq = wt_seq[:var_start] + mt_aa_change
                        if not mt_stop_codon:
                            mt_seq += wt_seq[var_start+len(wt_aa_change):]
                        

                        if (len(wt_aa_change) != 0 and len(mt_aa_change) != 0):
                            if wt_aa_change[0] == mt_aa_change[0]:
                                wt_aa_change = wt_aa_change[1:]
                                mt_aa_change = mt_aa_change[1:]
                                var_start += 1
                        

                    # print(f'----------{csq}---------')
                    # print(f'var_start: {var_start}')
                    # print(f'wt_seq: {wt_seq}')
                    # print(f'mt_seq: {mt_seq}')

                    variant_effects.change_entry(chrom=chrom,
                                               start=start,
                                               end=stop,
                                               gene_id=gene_id,
                                               gene_name=gene_name,
                                               transcript_id=transcript_id,
                                               transcript=None,
                                               transcript_bp=None,
                                               source=entry.INFO["SRC"],
                                               group=entry.INFO["GRP"],
                                               var_type=csq,
                                               var_start=var_start,
                                               wt_seq=wt_seq,
                                               mt_seq=mt_seq,
                                               vaf=vaf,
                                               ao=mt_ad,
                                               dp=dp,
                                               nmd_event=nmd)

                    # check if variant is self dissimilar
                    if variant_effects.self_dissimilarity():
                        variant_effects.write_entry()



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


    @staticmethod
    def get_variant_startpos(protein_position):
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
            consequence = "inframe_INS"
        elif "inframe_deletion" in consequences:
            consequence = "inframe_DEL"
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
        if (entry.INFO['SRC'] == "long_indel" or 
            entry.INFO['SRC'] == "exitron" or 
            entry.INFO['SRC'] == "alt_5prime" or 
            entry.INFO['SRC'] == "alt_3prime" or
            entry.INFO['SRC'] == "intron_retention" or 
            entry.INFO['SRC'] == "exon_skip" or 
            entry.INFO['SRC'] == "mutex_exons"):

            mt_ad = entry.INFO['AO']
            #wt_ad = str(int(float(entry.INFO['DP']) - float(mt_ad)))

        return mt_ad
                            

# augmented transcript
class VariantEffects:
    def __init__(self, 
                 output_dir, 
                 anno_file,
                 counts_file):
        self.anno = reference.Annotations(anno_file).anno
        self.counts = reference.Counts(counts_file).counts
        self.data = {}

        output_dir_path =  Path(output_dir)
        if not output_dir_path.exists():
            output_dir_path.mkdir()

        self.variantEffectsFile = Path(output_dir_path, "variant_effects.tsv")
        # create dir if it not exists

        self.fh = open(self.variantEffectsFile, 'w')


    def change_entry(self, chrom, start, end, gene_id, gene_name, 
                     transcript_id, transcript, transcript_bp, source, 
                     group, var_type, var_start, wt_seq, mt_seq, vaf, ao, 
                     dp, nmd_event):

        self.data = {}
        self.data["chrom"] = chrom
        self.data["start"] = start
        self.data["end"] = end
        self.data["gene_id"] = gene_id
        self.data["gene_name"] = gene_name
        self.data["transcript_id"] = transcript_id
        self.data["transcript"] = transcript
        self.data["transcript_bp"] = transcript_bp
        self.data["source"] = source
        self.data["group"] = group
        self.data["var_type"] = var_type
        self.data["var_start"] = var_start
        self.data["wt_seq"] = VariantEffects.adjust_wildtype(wt_seq, mt_seq)
        self.data["mt_seq"] = mt_seq
        self.data["vaf"] = vaf
        self.data["ao"] = ao
        self.data["dp"] = dp

        var_bnds = self.determine_var_bnds()
        self.data["aa_var_start"], self.data["aa_var_end"] = var_bnds

        self.determine_subsequence()
        self.determine_NMD(nmd_event)

        # determine TPM
        self.get_counts()



    @staticmethod
    def adjust_wildtype(wt, mt):
        """ make sure that wildtype sequence is of same length as mutant"""
        wt_seq = wt
        if len(mt) >= len(wt):
            for i in range(len(wt), len(mt)):
                wt_seq += '$'
        else:
            wt_seq = wt[:len(mt)]
        return wt_seq
    
    def determine_var_bnds(self):
        """determines the actual variant start and end by comparing
        the wildtype and mutant sequences"""
        variants = []
        for i in range(self.data["var_start"], len(self.data["mt_seq"])):
            if self.data["wt_seq"][i] != self.data["mt_seq"][i]:
                variants.append(i)

        if len(variants) != 0: # 
            return min(variants), max(variants)
        else:
            return -1, -1
    
    def determine_subsequence(self):
        """determines the subsequence of the variant"""

        shift = 0
        # determine start of subsequence
        if self.data["aa_var_start"] <= 24:
            left = 0
            new_var_start = self.data["aa_var_start"]
        else:
            left = self.data["aa_var_start"] - 24
            new_var_start = 24

        # determine end of subsequence
        if self.data["aa_var_end"] + 24 > len(self.data["mt_seq"]):
            right = len(self.data["mt_seq"]) - 1
        else:
            right = self.data["aa_var_end"] + 24
        var_len = self.data["aa_var_end"] - self.data["aa_var_start"] + 1
        new_var_end = new_var_start + var_len - 1

        """shift aa_var_[start|end] to the left to keep the correct position
        in the subsequence"""
        self.data["aa_var_start"] = new_var_start
        self.data["aa_var_end"] = new_var_end

        self.data["wt_subseq"] = self.data["wt_seq"][left:right+1]
        self.data["mt_subseq"] = self.data["mt_seq"][left:right+1]

    # TODO: check for NMD out of VEP
    def determine_NMD(self, nmd_event):
        if nmd_event is not None:
            self.data["NMD"] = nmd_event
        else:
            self.data["NMD"] = None

        self.data["PTC_exon_number"] = None
        self.data["PTC_dist_ejc"] = None
        self.data["NMD_escape_rule"] = None

        transcript = self.data["transcript"]
        # only need to be considered for framshift events
        if "frameshift" in self.data["var_type"]:
            # transcript is required to determine NMD
            if transcript is not None: 
                transcript_bp = int(self.data["transcript_bp"])
                start = self.data["start"]

                """ determine the coding sequence and adjust the breakpoint 
                (considering the adjusted startpos, e.g., start codon)"""
                cds, cds_bp = self.determine_cds(transcript,
                                                 transcript_bp)

                if cds is not None: # when the start codon could be found
                    """determine the genomic coordinate of the breakpoint -
                    consider start of segment 2 in fusion transcripts"""
                    """determine the genomic coordinate of the cds breakpoint 
                    - note: in fusion events this the start of the segment 2. 
                    We can use this since self.data["transcript_bp"] corresponds 
                    to seg2 in start"""
                    if self.data["source"] == "fusion":
                        bp_coord = int(start.split('|')[1])
                    else:
                        bp_coord = int(start) + transcript_bp + 1

                    # search for stop_codon
                    stop_pos, stop_coord = self.find_stop_codon(cds,
                                                                cds_bp,
                                                                bp_coord)

                    # check of stop_codon is PTC
                    if stop_pos != -1:
                        if self.data["transcript_id"] is not None:
                            if self.data["source"] == "fusion":
                                """ needs to be tid in 2nd transcript - mainly
                                because this is where the PTC is supposed to be"""
                                tid = self.data["transcript_id"].split('|')[1]
                                if tid == '.':
                                    return
                            else:
                                tid = self.data["transcript_id"]

                            exoninfo = self.anno[tid]
                            exons = list(exoninfo.keys())

                            exon_num, dist_ejc = self.annotate_stop_codon(exoninfo,
                                                                          stop_coord)

                            if exon_num != -1:
                                self.data["PTC_exon_number"] = f'{exon_num}'
                                if max(exons) > 1:
                                    self.data["PTC_dist_ejc"] = dist_ejc

                                nmd_escape = self.check_escape(exoninfo,
                                                               stop_coord,
                                                               exon_num,
                                                               dist_ejc)

                                if nmd_escape != -1:
                                    self.data["NMD"] = "NMD_escaping_variant"
                                    self.data["NMD_escape_rule"] = nmd_escape
                                else:
                                    self.data["NMD"] = "NMD_variant"


    
    @staticmethod
    def determine_cds(transcript, transcript_bp):
        """ determines the coding sequence of the transcript - 
        return coding sequence (cds) and its breakpoint within the cds
        (cds_bp) and its genomic coordinate (cds_bp_global)"""

        start_codon = re.search(r'ATG', transcript)
        if start_codon:
            """start codon has to occur before the breakpoint 
            - note the breakpoint corresponds to the local site"""
            if start_codon.start() < transcript_bp:
                # determines the cds - start after start codon...
                cds = transcript[start_codon.start()+3:]
                # ... and the breakpoint position within the cds
                cds_bp = transcript_bp - (start_codon.start() + 3)
                # """determine the genomic coordinate of the cds breakpoint 
                # - note: in fusion events this the start of the segment 2. 
                # We can use this since self.data["transcript_bp"] corresponds to 
                # start(seg1|seg2)"""
                # start = self.data["start"]
                # if self.data["source"] == "fusion":
                    # # in fusion events this contains the starts of both segments
                    # cds_bp_coord = start.split('|')[1]
                # else:
                    # cds_cp_coord = start + self.data["transcript_bp"] + 1
                return cds, cds_bp
            else:
                return None, None
        else:
            return None, None


    @staticmethod
    def find_stop_codon(cds, bp, bp_coord):
        """search for stop codon in coding sequence returns 0-based index

        the stop codon and 
        genomic coordinate of the end of the second segment
        """
        stop_pos = -1
        for i in range(0, len(cds), 3):
            codon = cds[i:i+3]
            if codon == "TAA" or codon == "TAG" or codon == "TGA":
                stop_pos = i # stop codon has been found
                break

        # stop codon shouldn't be found before the actual mutant in transcript
        if stop_pos <= bp:
            return -1, -1
        else:
            return stop_pos, bp_coord + (stop_pos - bp)  


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


    def self_dissimilarity(self):
        """ check if wildtype and mutant sequence are dissimilar """
        if self.data["wt_subseq"] != self.data["mt_subseq"]:
            return True
        else:
            return False

    def get_counts(self):
        """ determine the counts for the variant"""

        # fusion events contain gene_id/chrom of both segments
        if self.data["source"] == "fusion":
            gene_ids = self.data["gene_id"].split('|')
            chroms = self.data["chrom"].split('|')
            key1 = (gene_ids[0], chroms[0])
            key2 = (gene_ids[1], chroms[1])

            tpm = ''
            if key1 in self.counts:
                tpm = f'{self.counts[key1][self.data["group"]]}|'
            else:
                tpm = f'NA|'

            if key2 in self.counts:
                tpm += f'{self.counts[key2][self.data["group"]]}'
            else:
                tpm += f'NA'

            self.data["TPM"] = tpm

        else:
            key = (self.data["gene_id"], self.data["chrom"])
            if key in self.counts:
                self.data["TPM"] = self.counts[key][self.data["group"]]
            else:
                self.data["TPM"] = None



    # TODO: no hard coding of output
    def write_header(self):
        self.fh.write("chrom\tstart\tend\tgene_id\tgene_name\ttranscript_id\t")
        self.fh.write("source\tgroup\tvar_type\twt_subseq\tmt_subseq\t")
        self.fh.write("var_start\taa_var_start\taa_var_end\t")
        self.fh.write("vaf\tao\tdp\tTPM\tNMD\tPTC_dist_ejc\tPTC_exon_number\t")
        self.fh.write("NMD_escape_rule\n")
        
    def write_entry(self):
        self.fh.write(f'{format_output(self.data["chrom"])}\t')
        self.fh.write(f'{format_output(self.data["start"])}\t')
        self.fh.write(f'{format_output(self.data["end"])}\t')
        self.fh.write(f'{format_output(self.data["gene_id"])}\t')
        self.fh.write(f'{format_output(self.data["gene_name"])}\t')
        self.fh.write(f'{format_output(self.data["transcript_id"])}\t')
        self.fh.write(f'{format_output(self.data["source"])}\t')
        self.fh.write(f'{format_output(self.data["group"])}\t')
        self.fh.write(f'{format_output(self.data["var_type"])}\t')
        self.fh.write(f'{format_output(self.data["wt_subseq"])}\t')
        self.fh.write(f'{format_output(self.data["mt_subseq"])}\t')
        self.fh.write(f'{format_output(self.data["var_start"])}\t')
        self.fh.write(f'{format_output(self.data["aa_var_start"])}\t')
        self.fh.write(f'{format_output(self.data["aa_var_end"])}\t')
        self.fh.write(f'{format_output(self.data["vaf"])}\t')
        self.fh.write(f'{format_output(self.data["ao"])}\t')
        self.fh.write(f'{format_output(self.data["dp"])}\t')
        self.fh.write(f'{format_output(self.data["TPM"])}\t')
        self.fh.write(f'{format_output(self.data["NMD"])}\t')
        self.fh.write(f'{format_output(self.data["PTC_dist_ejc"])}\t')
        self.fh.write(f'{format_output(self.data["PTC_exon_number"])}\t')
        self.fh.write(f'{format_output(self.data["NMD_escape_rule"])}\n')
        
    def close_file(self):
        self.fh.close()

def format_output(field):
    if field == None:
        return '.'
    else:
        return field

# deter
def find_all_lowercase(seq):
    matches = []
    for i in range(len(seq)):
        if seq[i].islower():
            matches.append(i)
    return matches
