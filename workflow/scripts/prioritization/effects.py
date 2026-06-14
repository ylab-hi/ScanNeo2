"""Derive the protein-level effect of each variant.

Reconstructs the wildtype/mutant peptide subsequence around each variant, determines NMD status and
PTC location for frameshift events, attaches expression (TPM), and writes the per-variant-type
`*_variant_effects.tsv` consumed by the binding-affinity prediction step.
"""

# classes
import utility as ut
import reference

# standard
import re
from pathlib import Path

# augmented transcript
class VariantEffects:
    def __init__(self, options, vartype):

        self.annotation = reference.Annotation(options.reference, options.anno)
        self.exome = self.annotation.exome
        self.counts = reference.Counts(options.counts).counts
        self.data = {}


        output_dir_path =  Path(options.output_dir)
        if not output_dir_path.exists():
            output_dir_path.mkdir()

        self.variantEffectsFile = Path(output_dir_path, f"{vartype}_variant_effects.tsv")
        # create dir if it not exists


        self.fh = open(self.variantEffectsFile, 'w')
        self.write_header()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.fh.close()
        return False


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

    def determine_NMD(self, nmd_event):
        if nmd_event is not None:
            self.data["NMD"] = nmd_event
        else:
            self.data["NMD"] = None

        self.data["PTC_exon_number"] = None
        self.data["PTC_dist_ejc"] = None
        self.data["NMD_escape_rule"] = None

        transcript = self.data["transcript"]
        if "frameshift" not in self.data["var_type"] or transcript is None:
            return
        if self.data["transcript_id"] is None:
            return

        transcript_bp = int(self.data["transcript_bp"])

        # Locate the canonical start codon's mRNA position. For non-fusion
        # variants the GTF-derived start_codons table is authoritative; for
        # fusions (arriba's pre-spliced chimeric peptide has no canonical GTF
        # start codon) keep the leading-ATG heuristic on the fusion sequence.
        if self.data["source"] == "fusion":
            start_codon_offset = transcript.find('ATG')
            if start_codon_offset == -1:
                return
        else:
            canonical_pos = self.annotation.start_codons.get(
                self.data["transcript_id"])
            if canonical_pos is None:
                return
            start_codon_offset = self.annotation.genomic_to_mrna(
                self.data["transcript_id"], canonical_pos)
            if start_codon_offset is None:
                return

        cds, cds_bp = self.determine_cds(transcript, transcript_bp,
                                          start_codon_offset)
        if cds is None:
            return

        stop_pos = self.find_stop_codon(cds, cds_bp)
        if stop_pos == -1:
            return

        if self.data["source"] == "fusion":
            # Fusion: arriba reports segment-2's genomic breakpoint as a 1-based
            # position; shift to 0-based to match the exon table convention.
            # The transcript is already spliced by arriba so linear arithmetic
            # past the breakpoint is correct within segment 2.
            tid = self.data["transcript_id"].split('|')[1]
            if tid == '.':
                return
            seg2_bp_0based = int(self.data["start"].split('|')[1]) - 1
            stop_coord = seg2_bp_0based + (stop_pos - cds_bp)
            # Arriba doesn't expose strand directly; assume +.
            strand = '+'
        else:
            tid = self.data["transcript_id"]
            # CDS position → mRNA position → genomic coord (strand-aware).
            mrna_stop = stop_pos + start_codon_offset + 3
            stop_coord = self.annotation.mrna_to_genomic(tid, mrna_stop)
            if stop_coord is None:
                return
            strand = self.annotation.transcriptome[tid][2]

        exoninfo = self.exome[tid]
        exons = list(exoninfo.keys())

        exon_num, dist_ejc = self.annotate_stop_codon(exoninfo, stop_coord, strand)
        if exon_num == -1:
            return

        self.data["PTC_exon_number"] = f'{exon_num}'
        if max(exons) > 1:
            self.data["PTC_dist_ejc"] = dist_ejc

        nmd_escape = self.check_escape(exoninfo, stop_coord, exon_num,
                                       dist_ejc, strand)
        if nmd_escape != -1:
            self.data["NMD"] = "NMD_escaping_variant"
            self.data["NMD_escape_rule"] = nmd_escape
        else:
            self.data["NMD"] = "NMD_variant"


    
    @staticmethod
    def determine_cds(transcript, transcript_bp, start_codon_offset):
        """Carve the CDS from the spliced mRNA given the known start codon
        position. The caller (`determine_NMD`) computes `start_codon_offset`:
        for non-fusion variants from the GTF-derived `Annotation.start_codons`
        table (the canonical CDS start, not the first ATG); for fusions from
        `transcript.find('ATG')` on arriba's already-spliced fusion peptide.

        Returns (cds, cds_bp). Returns (None, None) when the start codon is
        absent or sits at or beyond the variant breakpoint.
        """
        if start_codon_offset is None or start_codon_offset >= transcript_bp:
            return None, None
        cds = transcript[start_codon_offset+3:]
        cds_bp = transcript_bp - (start_codon_offset + 3)
        return cds, cds_bp


    @staticmethod
    def find_stop_codon(cds, bp):
        """Return the 0-based CDS-relative position of the first in-frame stop
        codon strictly downstream of bp, or -1 if none. The caller is
        responsible for converting this CDS position into a genomic coordinate
        (strand-aware for non-fusion; linear arithmetic for fusion)."""
        stop_pos = -1
        for i in range(0, len(cds), 3):
            codon = cds[i:i+3]
            if codon == "TAA" or codon == "TAG" or codon == "TGA":
                stop_pos = i
                break

        if stop_pos <= bp:
            return -1
        return stop_pos


    def annotate_stop_codon(self, exoninfo, stop_coord, strand):
        """Locate which exon holds the PTC and report `dist_ejc`, the number of
        bases between the PTC and the 3' end of that exon in transcript
        orientation. Coordinates are 0-based half-open."""
        for exon_number in exoninfo:
            exon_start, exon_end = exoninfo[exon_number]
            if exon_start <= stop_coord < exon_end:
                if strand == '+':
                    dist_ejc = (exon_end - 1) - stop_coord
                else:
                    dist_ejc = stop_coord - exon_start
                return exon_number, dist_ejc

        return -1, -1


    def check_escape(self, exoninfo, stop_coord, exon_num, dist_ejc, strand):
        """Determine whether a PTC escapes NMD. Returns one of:

          1 — Intronless transcript (only one exon). No downstream
              exon-exon junction → no EJC for the NMD machinery to recruit.
          2 — Start-proximal PTC, within 100 nt of the 5' end of exon 1 in
              transcript orientation. Translation re-initiation downstream
              can rescue the transcript.
          3 — PTC within 50 nt of the last exon-exon junction (i.e. close to
              the 3' end of the penultimate exon in transcript orientation).
              The downstream EJC is removed during the first round of
              translation before NMD can act.
          4 — PTC in the last exon — no downstream EJC at all.
         -1 — None of the above; NMD targets the transcript.

        ASCII layout (transcript orientation, 5'→3'):

          rule 2:   vvv
                ..ES...EE..I.ES...EE.I.ES....EE.I.ES....EE
          rule 3:                       vvv
                  ES...EE..I.ES...EE.I.ES....EE.I.ES....EE
          rule 4:                                    vvvv
                  ES...EE..I.ES...EE.I.ES....EE.I.ES....EE

        (ES = exon 5' end, EE = exon 3' end, I = intron, v = PTC.)
        """
        exons = list(exoninfo.keys())
        last_exon = int(max(exons))

        if exon_num == 1:
            if last_exon == 1: # there is only one exon
                return 1
            elif last_exon > 1:
                exon_start, exon_end = exoninfo[1]
                # rule 2: PTC within 100 nt of the 5' end of exon 1 in
                # transcript orientation
                if strand == '+':
                    dist_from_5p = stop_coord - exon_start
                else:
                    dist_from_5p = (exon_end - 1) - stop_coord
                if dist_from_5p < 100:
                    return 2
        if exon_num == last_exon-1:
            # rule 3: PTC within 50 nt of the 3' end of the penultimate exon
            # in transcript orientation — i.e. close to the last EJC
            if dist_ejc <= 50:
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
                if self.data["group"] in self.counts[key]:
                    self.data["TPM"] = self.counts[key][self.data["group"]]
                else:
                    self.data["TPM"] = None
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
