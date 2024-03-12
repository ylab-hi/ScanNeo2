# classes
import utility as ut
import effects
import reference

# standard
import re
from pathlib import Path

class Fusions:
    def __init__(self, 
                 fusions_input,
                 options,
                 vartype):

        proteome = reference.Proteome(options.proteome).proteome
        annotation = reference.Annotation(options.reference, options.anno)

        # create variant effects object
        self.variant_effects = effects.VariantEffects(options, vartype)

        group = Path(fusions_input).stem.split('_fusions')[0]
        confidence_level = {"low": ["low", "medium", "high"],
                           "medium": ["medium", "high"],
                           "high": ["high"]}


        with open(fusions_input) as fh:
            next(fh) # ignore header
            for line in fh:
                cols = line.rstrip().split('\t')
                if cols[14] in confidence_level[options.confidence]:
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

                    self.variant_effects.change_entry(chrom=chrom,
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
                    if self.variant_effects.self_dissimilarity():
                        self.variant_effects.write_entry()


            self.variant_effects.close_file()


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
            variants = ut.find_all_lowercase(fusion_seg1)
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
