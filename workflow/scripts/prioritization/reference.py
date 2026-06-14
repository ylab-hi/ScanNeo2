"""Load reference data used throughout prioritization.

`Proteome` indexes the reference peptide FASTA by transcript ID, `Annotation` parses the GTF into
transcriptome and exome coordinate maps, and `Counts` loads the per-gene expression count table.
"""

import pyfaidx
import re

COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]

class Proteome:
    def __init__(self, fastaFile):
        # parse proteome
        ref = pyfaidx.Fasta(fastaFile)
        self.proteome = {}
        for record in ref:
            identifier = record.long_name
            match = re.search(r'transcript:([^.\s]+)', identifier)
            if match:
                tid = match.group(1)
                self.proteome[tid] = str(record)



class Annotation:
    def __init__(self, fastaFile, gtfFile):
        self.ref = pyfaidx.Fasta(fastaFile)
        self.transcriptome = self.parse_transcriptome(gtfFile)
        self.exome = self.parse_exome(gtfFile)
        self.exon_map_cache = {}

    # GTF coords are 1-based inclusive; we store 0-based half-open throughout
    # (start inclusive, end exclusive) for direct compatibility with pyfaidx
    # slicing and Python range semantics.
    def parse_transcriptome(self, gtfFile):
        transcriptome = {}
        with open(gtfFile, "r") as fh:
            for line in fh:
                # ignore header
                if not line.startswith("#"):
                    l = line.rstrip().split("\t")
                    if l[2] == "transcript":
                        transcript_id = re.search(r'transcript_id "([^.\s]+)"', l[8]).group(1)
                        if transcript_id not in transcriptome:
                            transcriptome[transcript_id] = [int(l[3])-1, int(l[4]), l[6]]

        return transcriptome


    def parse_exome(self, gtfFile):
        exome = {}
        with open(gtfFile, 'r') as fh:
            for line in fh:
                # ignore header
                if not line.startswith('#'):
                    l = line.rstrip().split('\t')

                    if l[2] == 'exon':
                        exon_number = re.search(r'exon_number "([^.\s]+)"', l[8]).group(1)
                        transcript_id = re.search(r'transcript_id "([^.\s]+)"', l[8]).group(1)

                        if transcript_id not in exome:
                            exome[transcript_id] = {}

                        if not exon_number in exome[transcript_id]:
                            exome[transcript_id][int(exon_number)] = [int(l[3])-1, int(l[4])]

        return exome


    def exon_map(self, transcript_id):
        """Per-exon table in 5'→3' (mRNA) order.

        GTF exon_number is by transcript position — exon 1 is always the first
        exon transcribed, regardless of strand. Iterating ascending therefore
        walks the mRNA from its 5' to 3' end on either strand.

        Returns (strand, exons) where exons is a list of tuples
        (mrna_start, mrna_end_excl, genomic_start, genomic_end_excl),
        all 0-based half-open.
        """
        if transcript_id in self.exon_map_cache:
            return self.exon_map_cache[transcript_id]

        strand = self.transcriptome[transcript_id][2]
        exoninfo = self.exome[transcript_id]
        exons = []
        mrna_cursor = 0
        for n in sorted(exoninfo.keys()):
            g_start, g_end = exoninfo[n]
            length = g_end - g_start
            exons.append((mrna_cursor, mrna_cursor + length, g_start, g_end))
            mrna_cursor += length

        result = (strand, exons)
        self.exon_map_cache[transcript_id] = result
        return result


    def splice_transcript(self, transcript_id, chrom):
        """Build the strand-aware spliced mRNA for a transcript.

        For - strand transcripts each exon's forward-strand genomic sequence is
        reverse-complemented before concatenation, so the returned mRNA reads
        5'→3' in transcript orientation.

        Returns (mrna, strand, exons) — exons as from exon_map.
        """
        strand, exons = self.exon_map(transcript_id)
        pieces = []
        for (_, _, g_start, g_end) in exons:
            seq = str(self.ref[chrom][g_start:g_end])
            if strand == '-':
                seq = revcomp(seq)
            pieces.append(seq)
        return ''.join(pieces), strand, exons


    def genomic_to_mrna(self, transcript_id, genomic_pos):
        """Map a 0-based genomic position to a 0-based mRNA position.

        Returns None if the position falls in an intron or outside the
        transcript. For - strand, the mRNA position counts from the 5' end of
        the transcript, which is the genomically-rightmost exon — so the
        genomically-rightmost base of an exon is the 5'-most base of that
        exon's mRNA chunk.
        """
        strand, exons = self.exon_map(transcript_id)
        for (m_start, _, g_start, g_end) in exons:
            if g_start <= genomic_pos < g_end:
                if strand == '+':
                    return m_start + (genomic_pos - g_start)
                else:
                    return m_start + (g_end - 1 - genomic_pos)
        return None


    def mrna_to_genomic(self, transcript_id, mrna_pos):
        """Inverse of genomic_to_mrna. Returns a 0-based genomic position, or
        None if mrna_pos is outside the transcript's spliced length."""
        strand, exons = self.exon_map(transcript_id)
        for (m_start, m_end, g_start, g_end) in exons:
            if m_start <= mrna_pos < m_end:
                offset = mrna_pos - m_start
                if strand == '+':
                    return g_start + offset
                else:
                    return g_end - 1 - offset
        return None


    def splice_with_variant(self, transcript_id, chrom,
                             affected_start, affected_end, allele):
        """Build the mutant spliced mRNA for a non-fusion variant.

        Coordinates are 0-based half-open (vcfpy `affected_start`/`affected_end`
        convention). `allele` is VEP's forward-strand `Allele` field: `'-'` for
        pure deletion, the inserted bases (without anchor) for insertion, the
        alt base(s) for SNV/MNV.

        Returns (mt_mrna, var_mrna_start) — `var_mrna_start` is the 0-based
        mRNA position of the first inserted/altered base in the mutant mRNA.
        Returns (None, None) if the variant has no clean exonic representation
        (intronic, or spans an exon-intron boundary).
        """
        mrna, strand, _ = self.splice_transcript(transcript_id, chrom)

        is_insertion = (allele != '-' and affected_start == affected_end)
        is_deletion  = (allele == '-')

        if is_insertion:
            # vcfpy reports affected_start = POS for insertions, i.e. the
            # position one-past the anchor base on + strand. The anchor itself
            # sits at affected_start - 1 (0-based genomic) and is unchanged.
            m_anchor = self.genomic_to_mrna(transcript_id, affected_start - 1)
            if m_anchor is None:
                return None, None
            alt = revcomp(allele) if strand == '-' else allele
            # On + strand the inserted bases follow the anchor in mRNA order;
            # on - strand they precede it (anchor's genomic neighbour at
            # affected_start lies at mRNA position m_anchor - 1).
            insert_at = m_anchor + 1 if strand == '+' else m_anchor
            return mrna[:insert_at] + alt + mrna[insert_at:], insert_at

        if is_deletion:
            # Anchor at affected_start is unchanged; deleted bases span
            # [affected_start + 1, affected_end).
            first_changed = affected_start + 1
            last_changed = affected_end - 1
        else:
            # SNV/MNV: replace [affected_start, affected_end) entirely.
            first_changed = affected_start
            last_changed = affected_end - 1

        m_a = self.genomic_to_mrna(transcript_id, first_changed)
        m_b = self.genomic_to_mrna(transcript_id, last_changed)
        if m_a is None or m_b is None:
            return None, None

        lo, hi = min(m_a, m_b), max(m_a, m_b)
        alt = '' if is_deletion else (revcomp(allele) if strand == '-' else allele)
        return mrna[:lo] + alt + mrna[hi + 1:], lo

class Counts:
    def __init__(self, countFile):
        # parse counts
        self.counts = {}
        if countFile is not None and countFile != "":
            with open(countFile, 'r') as count_fh:
                lines = count_fh.readlines()
            groups = lines[0].rstrip().split('\t')[3:]
            for line in lines[1:]:
                cols = line.rstrip().split('\t')

                gene_id = cols[0]
                chrom = cols[1]

                key = (gene_id, chrom)

                if key not in self.counts:
                    self.counts[key] = {}
                    for group in groups:
                        self.counts[key][group] = float(cols[3+groups.index(group)])

