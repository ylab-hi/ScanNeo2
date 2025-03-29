# classes
import utility as ut
import reference

# standard
import re
from pathlib import Path

# augmented transcript
class VariantEffects:
    def __init__(self, options, vartype):

        self.exome = reference.Annotation(options.reference, options.anno).exome
        self.counts = reference.Counts(options.counts).counts
        self.data = {}


        output_dir_path =  Path(options.output_dir)
        if not output_dir_path.exists():
            output_dir_path.mkdir()

        self.variantEffectsFile = Path(output_dir_path, f"{vartype}_variant_effects.tsv")
        # create dir if it not exists

        self.fh = open(self.variantEffectsFile, 'w')
        self.write_header()


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

                            exoninfo = self.exome[tid]
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
        self.fh.write(f'{ut.format_output(self.data["chrom"])}\t')
        self.fh.write(f'{ut.format_output(self.data["start"])}\t')
        self.fh.write(f'{ut.format_output(self.data["end"])}\t')
        self.fh.write(f'{ut.format_output(self.data["gene_id"])}\t')
        self.fh.write(f'{ut.format_output(self.data["gene_name"])}\t')
        self.fh.write(f'{ut.format_output(self.data["transcript_id"])}\t')
        self.fh.write(f'{ut.format_output(self.data["source"])}\t')
        self.fh.write(f'{ut.format_output(self.data["group"])}\t')
        self.fh.write(f'{ut.format_output(self.data["var_type"])}\t')
        self.fh.write(f'{ut.format_output(self.data["wt_subseq"])}\t')
        self.fh.write(f'{ut.format_output(self.data["mt_subseq"])}\t')
        self.fh.write(f'{ut.format_output(self.data["var_start"])}\t')
        self.fh.write(f'{ut.format_output(self.data["aa_var_start"])}\t')
        self.fh.write(f'{ut.format_output(self.data["aa_var_end"])}\t')
        self.fh.write(f'{ut.format_output(self.data["vaf"])}\t')
        self.fh.write(f'{ut.format_output(self.data["ao"])}\t')
        self.fh.write(f'{ut.format_output(self.data["dp"])}\t')
        self.fh.write(f'{ut.format_output(self.data["TPM"])}\t')
        self.fh.write(f'{ut.format_output(self.data["NMD"])}\t')
        self.fh.write(f'{ut.format_output(self.data["PTC_dist_ejc"])}\t')
        self.fh.write(f'{ut.format_output(self.data["PTC_exon_number"])}\t')
        self.fh.write(f'{ut.format_output(self.data["NMD_escape_rule"])}\n')
        
    def close_file(self):
        self.fh.close()
